#include "mmbam/mpileup.h"
#include "mmbam/mbgzf.h"
#include "mmbam/bam.h"
#include "mmbam/index.h"
#include <string>
#include <list>
#include <functional>

constexpr size_t max_buff = 256 * 1024 * 1024;

const bam_rec_t* multi_buffer_next_read(
        std::vector<std::vector<uint8_t>>& m_buffer,
        std::vector<size_t>& m_buffer_offsets,
        size_t& buffer_id) {
    size_t which_buffer = 0;
    int32_t smallest_pos = 0;
    const bam_rec_t* first_read = nullptr;
    for(size_t i_buffer = 0; i_buffer<m_buffer.size(); i_buffer++) {
        auto offset = m_buffer_offsets[i_buffer];
        if(offset < m_buffer[i_buffer].size()) {
            // i_buffer has read at offset
            const bam_rec_t *read = BAMREF(m_buffer[i_buffer].data() + offset);
            if(first_read == nullptr || read->pos < smallest_pos) {
                which_buffer = i_buffer;
                smallest_pos = read->pos;
                first_read = read;
            }
        }
    }

    if(first_read != nullptr) {
        // read found, increment buffer offset
        m_buffer_offsets[which_buffer] += (first_read->block_size + 4);
    }

    buffer_id = which_buffer;


    return first_read;
}

mpileup_t* allocate_mpileups(const std::vector<std::vector<uint8_t>>& m_buffers) {
    mpileup_t *pileups = new mpileup_t[400];
    for(size_t i=0; i<400; i++) pileups[i].allocate(m_buffers);
    return pileups;
}

void release_mpileups(mpileup_t* pileups) {
    for(size_t i=0; i<400; i++) pileups[i].release();
    delete [] pileups;
}

void mpileup(std::vector<std::reference_wrapper<const mfile_t::ptr_t>> mfiles,
        std::vector<std::reference_wrapper<const index_t>> indices,
        uint32_t ref_id, int32_t pos_start, int32_t pos_end,
        const std::function<bool(const bam_rec_t&)>& predicate_func,
        const std::function<bool(const mpileup_t&)>& visitor_func) {

    const auto n_files = mfiles.size();

    if(n_files != indices.size()) throw std::runtime_error("number of files and indices mismatch");

    // load reads in region from all files
    std::vector<std::vector<uint8_t>> m_buffers(n_files);
    std::vector<size_t> m_buffer_offsets(n_files, 0);
    for(size_t i=0; i<n_files; i++) {
        m_buffers[i] = bam_load_region(mfiles[i], indices[i], ref_id, pos_start, pos_end);
    }
    bool debug = false;

    if(pos_start <= 15635 && pos_end > 15635) {
        debug = true;
    }

    // initialize pileup data structure
    auto pileups = allocate_mpileups(m_buffers);
    int p_head = 0;
    int p_tail = 0;
    int p_size = 400;

    // pileup until reads are exhausted
    size_t buffer_id;
    const bam_rec_t *bam_it;
    while((bam_it = multi_buffer_next_read(m_buffers, m_buffer_offsets, buffer_id)) != nullptr) {
        if(!READ_IN_REGION(bam_it, pos_start, pos_end)) { continue; }
        if(!predicate_func(*bam_it)) { continue; }

        while(p_head != p_tail) {
            if(pileups[p_head].ref_id < bam_it->ref_id || pileups[p_head].pos < bam_it->pos) {
                if(!visitor_func(pileups[p_head])) {
                    release_mpileups(pileups);
                    return;
                }
                p_head++;
                if(p_head == p_size) p_head = 0;
            }
            else break;
        }

        int p_iter = p_head;

        // walk through the read, and update pileups
        const uint32_t *cigar_ops = (const uint32_t *)(
                (const uint8_t *)(bam_it)
                + sizeof(bam_rec_t) 
                + sizeof(char) * bam_it->l_read_name );

        const uint8_t *bam_begin = m_buffers[buffer_id].data();
        size_t read_offset = BYTEREF(bam_it) - bam_begin;
        uint32_t qpos = 0;
        int32_t rpos = bam_it->pos;

        for(int op=0; op < bam_it->n_cigar_op; op++) {
            auto cigar_len = cigar_ops[op] >> 4;
            auto cigar_op  = cigar_ops[op] & 0xF;
            switch(cigar_op) {
                case 1: // I: insertion only consumes query 
                case 4: // S: soft padded only consumes query
                    qpos += cigar_len;
                    break;
                case 2: // D: deletion only consumes reference
                case 3: // N: reference skip only consumes reference
                    while(cigar_len > 0) {
        
                        if(p_iter != p_tail) {
                            mpileup_t *p = &pileups[p_iter];
                            if(p->depth[buffer_id] < PILEUP_MAX_DEPTH) p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                .buffer_id = buffer_id,
                                .offset = read_offset, 
                                .qpos = qpos,
                                .is_deletion = true};
                            p_iter++;
                            if(p_iter == p_size) p_iter = 0;
                        }
                        else {
                            int next_tail = p_tail + 1;
                            if(next_tail == p_size) next_tail = 0;
                            if(next_tail != p_head) {
                                mpileup_t *p = &pileups[p_tail];
                                p->set_location(bam_it->ref_id, rpos, m_buffers);
                                p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                    .buffer_id = buffer_id,
                                    .offset = read_offset,
                                    .qpos = qpos,
                                    .is_deletion = true};
                                p_tail = next_tail;
                                p_iter = p_tail;
                            } else {
                                std::cerr<<"circular buffer full, result truncated"<<std::endl;
                            }
                        }

                        rpos++;
                        cigar_len--;
                    }
                    break;
                case 0: // M: alignment match consumes both
                case 7: // =: sequence match consumes both
                case 8: // X: sequence mismatch consumes both

                    while(cigar_len > 0) {
        
                        if(p_iter != p_tail) {
                            mpileup_t *p = &pileups[p_iter];
                            //if(p->depth < PILEUP_MAX_DEPTH) p->info[p->depth++] = {read_offset, qpos, false};
                            if(p->depth[buffer_id] < PILEUP_MAX_DEPTH) {
                                p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                .buffer_id = buffer_id,
                                .offset = read_offset, 
                                .qpos = qpos,
                                .is_deletion = false};
                            }
                            p_iter++;
                            if(p_iter == p_size) p_iter = 0;
                        }
                        else {
                            int next_tail = p_tail + 1;
                            if(next_tail == p_size) next_tail = 0;
                            if(next_tail != p_head) {
                                mpileup_t *p = &pileups[p_tail];
                                p->set_location(bam_it->ref_id, rpos, m_buffers);
                                p->get_info(buffer_id, p->depth[buffer_id]++) = {
                                    .buffer_id = buffer_id,
                                    .offset = read_offset,
                                    .qpos = qpos,
                                    .is_deletion = false};
                                p_tail = next_tail;
                                p_iter = p_tail;
                            }
                            else {
                                std::cerr<<"circular buffer full, result truncated"<<std::endl;
                            }
                        }

                        qpos++;
                        rpos++;
                        cigar_len--;
                    }
                    break;
                default: // other operations consume neither
                    break;
            }
        }
    }

    while(p_head != p_tail) {
        // emit the rest
        if(p_head == p_size) p_head = 0;
        if(p_head == p_tail) break;

        if(pileups[p_head].pos >= pos_end) {
            release_mpileups(pileups);
            return;
        }

        if(!visitor_func(pileups[p_head++])) {
            release_mpileups(pileups);
            return;
        }
    }

    release_mpileups(pileups);
}
