#include "mmbam/mbgzf.h"
#include "mmbam/bam.h"
#include "mmbam/index.h"
#include <vector>
#include <iostream>
#include <cassert>

#include <tbb/parallel_for.h>

#define CHUNK_SIZE 16384

bool bam_buffer_contains_header(const std::vector<uint8_t>& buffer) {
    // fixed part of the header
    size_t curr_size = 8;
    if(buffer.size() < curr_size) return false;
    const uint8_t *buffptr = &buffer[0];

    // text and n_ref
    uint32_t l_text = *(uint32_t*)(buffptr+4);
    curr_size += (l_text + 4);
    if(buffer.size() < curr_size) return false;

    uint32_t n_ref = *(uint32_t*)(buffptr + l_text + 8);
    for(uint32_t i_ref = 0; i_ref<n_ref; ++i_ref) {
        // need l_name
        curr_size += 4;
        if(buffer.size() < curr_size) return false;
        uint32_t l_name = *(uint32_t *)(buffptr + curr_size - 4);
        // need name and l_ref
        curr_size += l_name + 4;
        if(buffer.size() < curr_size) return false;
    }
    return true;
};

std::vector<uint8_t> bam_load_block(const mfile_t::ptr_t& mfile, uint64_t ioffset_first, uint64_t ioffset_last) {

    if(ioffset_first == ioffset_last) return std::vector<uint8_t>();

    auto coffset_first = index_coffset(ioffset_first);
    auto uoffset_first = index_uoffset(ioffset_first);

    auto coffset_last = index_coffset(ioffset_last);
    auto uoffset_last = index_uoffset(ioffset_last);

    auto bgzf_head        = bgzf_iterator_at(mfile, 0);
    auto bgzf_block_first = bgzf_iterator_at(mfile, coffset_first);
    auto bgzf_block_last  = bgzf_iterator_at(mfile, coffset_last);

    // special case: coffsets are the same
    if(coffset_first == coffset_last) {
        auto inflated_bytes = bgzf_inflate(*bgzf_block_first);

        auto first_byte = inflated_bytes.cbegin() + uoffset_first;
        auto last_byte  = uoffset_last == uoffset_first ? 
            inflated_bytes.cend() :
            inflated_bytes.cbegin() + uoffset_last;

        return std::vector<uint8_t>(first_byte, last_byte);
    }


    // inflate from bgzf_block_first to bgzf_block_last
    // trim uoffset_first from the beginning
    // trim uoffset_last from the end

    auto bgzf_block_it = bgzf_block_first;

    //std::vector<uint8_t> inflated_bytes;
    //inflated_bytes = bgzf_inflate_range(
    //        begin<const uint8_t>(mfile) + coffset_first,
    //        coffset_last - coffset_first);

    auto inflated_bytes = bgzf_inflate_range_p(
            begin<const uint8_t>(mfile) + coffset_first,
            coffset_last - coffset_first, [](
        auto* src, auto src_len, auto dest_len,
        auto& src_off_vector, auto& dest_off_vector,
        auto inflate) -> std::vector<uint8_t> {

        std::vector<uint8_t> buffer;
        buffer.resize(dest_len);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, src_off_vector.size()), [&](const auto& r){

            auto sl = ( r.end() == src_off_vector.size() ? src_len : src_off_vector[r.end()] ) - src_off_vector[r.begin()];
            auto dl = ( r.end() == dest_off_vector.size() ? dest_len : dest_off_vector[r.end()] ) - dest_off_vector[r.begin()];

            inflate(src+src_off_vector[r.begin()], sl, &buffer[dest_off_vector[r.begin()]], dl);

        });

        return buffer;
    });





    if(coffset_last < mfile->size) {
        auto last_block = bgzf_inflate(*bgzf_block_last);
        inflated_bytes.insert(inflated_bytes.end(), last_block.cbegin(), last_block.cbegin() + uoffset_last);
    }

    return std::vector<uint8_t>(inflated_bytes.cbegin() + uoffset_first, inflated_bytes.cend());
    
}

size_t bam_count_records(const std::vector<uint8_t>& buffer) {
    auto bam_it = bam_iterator(reinterpret_cast<const bam_rec_t*>(&buffer[0]));
    auto bam_it_beg = bam_iterator(reinterpret_cast<const bam_rec_t*>(&buffer[0]));
    auto bam_it_end = bam_iterator(reinterpret_cast<const bam_rec_t*>(&buffer[0] + buffer.size()));
    size_t total_recs = 0;
    while(bam_it < bam_it_end) {
        total_recs++;
        bam_it++;
    }

    return total_recs;
}

int16_t bam_query_length(const bam_rec_t* b) {
	uint32_t ql = 0;
	const uint32_t *cigar_ops = (const uint32_t *)(
			(const uint8_t *)(b)
			+ sizeof(bam_rec_t)
			+ sizeof(char) * b->l_read_name );
	for(int op = 0; op < b->n_cigar_op; op++) {
		switch(cigar_ops[op] & 0xF) {
			case 0: // M
			case 2: // D
			case 7: // =
			case 8: // X
				ql += cigar_ops[op] >> 4;
				break;
			default:
				break;
		}
	}
	return ql;
}

std::vector<uint8_t> bam_load_region(const mfile_t::ptr_t& mfile, const index_t& index, int32_t ref_id, int32_t region_start, int32_t region_end) {

    if(ref_id >= index.n_ref) throw std::runtime_error("reference not found");
    if(index.ref[ref_id].n_intv == 0) throw std::runtime_error("empty reference");

    using intv_idx_t = decltype(index.ref[ref_id].n_intv);

    // start_interval
    intv_idx_t start_intv = region_start / CHUNK_SIZE;
    intv_idx_t end_intv   = region_end   / CHUNK_SIZE;

    if(start_intv > index.ref[ref_id].n_intv - 1)
        start_intv = index.ref[ref_id].n_intv - 1;

    while(end_intv < index.ref[ref_id].n_intv && 
            index.ref[ref_id].ioffset[start_intv] == index.ref[ref_id].ioffset[end_intv])
        end_intv++;

    if(end_intv > index.ref[ref_id].n_intv - 1)
        end_intv = index.ref[ref_id].n_intv - 1;

    // load bgzf block in batch
    std::vector<uint8_t> bam_buffer_preload = bam_load_block(
            mfile,
            index.ref[ref_id].ioffset[start_intv],
            index.ref[ref_id].ioffset[end_intv]);

    auto bam_iter = reinterpret_cast<const bam_rec_t *>(bam_buffer_preload.data());
    auto buffer_end = bam_buffer_preload.data() + bam_buffer_preload.size();
    size_t first_in_region = -1;

    if(bam_buffer_preload.size() > 0) {
        auto bam_next = BYTEREF(bam_iter) + bam_iter->block_size + 4;
        // reads returned. See if the region is contained
        while(bam_next < buffer_end) {
            if( READ_IN_REGION(bam_iter, region_start, region_end) && first_in_region == -1) {
                first_in_region = BYTEREF(bam_iter) - bam_buffer_preload.data();
            }
            if(bam_iter->ref_id != ref_id || bam_iter->pos >= region_end) {
                if(first_in_region == -1) return std::vector<uint8_t>();
                return std::vector<uint8_t>(
                        bam_buffer_preload.cbegin() + first_in_region,
                        bam_buffer_preload.cbegin() + (BYTEREF(bam_iter) - bam_buffer_preload.data())
                        );
            }
            bam_iter = BAMREF(bam_next);
            bam_next = BYTEREF(bam_iter) + bam_iter->block_size + 4;
        }
    }

    // last read is still within [region_start, region_end)
    // append data starting from end_intv
    auto ioffset = index.ref[ref_id].ioffset[end_intv];
    auto bgzf_it = bgzf_iterator_at(mfile, index_coffset(ioffset));
    auto bgzf_end = bgzf_iterator_at(mfile, mfile->size);

    std::vector<uint8_t> bam_buffer;

    bool found_outbound = false;

    size_t prev_offset = index_uoffset(ioffset);

    while(bgzf_it < bgzf_end && !found_outbound) {
        auto inflated = bgzf_inflate(*bgzf_it);
        bgzf_it++;

        bam_buffer.insert(bam_buffer.end(), inflated.cbegin(), inflated.cend());

        bam_iter = reinterpret_cast<const bam_rec_t *>(bam_buffer.data() + prev_offset);
        decltype(bam_iter) bam_prev;
        auto buffer_end = bam_buffer.data() + bam_buffer.size();
        while(BYTEREF(bam_iter) < buffer_end) {
            if( READ_IN_REGION(bam_iter, region_start, region_end) && first_in_region == -1 ) {
                first_in_region = BYTEREF(bam_iter) - bam_buffer.data() - index_uoffset(ioffset);
            }
            if(bam_iter->ref_id != ref_id || bam_iter->pos >= region_end) {
                found_outbound = true;
                break;
            }
            bam_prev = bam_iter;
            bam_iter = BAM_NEXT(bam_iter);
        }

        // did not find out-of-bound read, readjust prev_offset
        prev_offset = BYTEREF(bam_prev) - bam_buffer.data();
    }
    
    // concatenate preload buffer with bam_buffer(offset with first uoffset
    bam_buffer_preload.insert(bam_buffer_preload.end(),
            bam_buffer.cbegin() + index_uoffset(ioffset),
            bam_buffer.cbegin() + (BYTEREF(bam_iter) - bam_buffer.data()));

    if(first_in_region == -1)
        return std::vector<uint8_t>();

    return std::vector<uint8_t>(
            bam_buffer_preload.cbegin() + first_in_region,
            bam_buffer_preload.cend());

}
