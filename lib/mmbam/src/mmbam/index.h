#ifndef NVPL_INDEX_H
#define NVPL_INDEX_H

#include <vector>
#include "mfile.h"

constexpr auto idx_header_size = sizeof(char[4]) + sizeof(uint32_t);

//! A pair of ioffset that defines a genomic region
using region = std::pair<uint64_t, uint64_t>;

//! this struct mirrors the chunk level structure of the BAM index format
struct index_chunk_t {
    uint64_t beg;
    uint64_t end;
};

//! this struct mirrors the bin level structure of the BAM index format
struct index_bin_t {
    uint32_t bin_id;
    uint32_t n_chunk;
    index_chunk_t *chunk;
};

//! this struct mirrors the reference level structure of the BAM index format
struct index_ref_t {
    uint32_t n_bin;
    index_bin_t *bin;
    uint32_t n_intv;
    uint64_t *ioffset;
};

//! this struct mirrors the top level structure of the BAM index format
struct index_t {
    char magic[4];
    uint32_t n_ref;
    index_ref_t *ref;
    uint64_t n_no_coor;
};

constexpr auto MAX_BIN = (((1<<18)-1)/7);

//! calculates the coffset from the given virtual offset
inline uint64_t index_coffset(uint64_t offset) { return offset >> 16; }

//! calculates the uoffset from the given virtual offset
inline uint16_t index_uoffset(uint64_t offset) { return static_cast<uint16_t>(offset); }

//! release the memory used by the internal structure of a loaded BAI file
void index_free(index_t& index);

//! Parse a BAI file stream, and create the corresponding index_t object
//! \param ins An input stream that contains the BAI file content
//! \return Fully initialized index_t representing the BAI file
template<class InStream>
index_t index_read(InStream&& ins) {
    index_t index;
    ins.read((char *)&index, idx_header_size);

    // allocate ref blocks
    index.ref = new index_ref_t[index.n_ref];

    // read ref blocks
    for(size_t ref_i = 0; ref_i < index.n_ref; ref_i++) {
        index_ref_t *ref = index.ref + ref_i;

        ins.read((char *)ref, sizeof(uint32_t));
        
        // allocate bins
        ref->bin = new index_bin_t[ref->n_bin];

        // read bins
        for(size_t bin_i = 0; bin_i < ref->n_bin; bin_i++) {
            index_bin_t *bin = ref->bin + bin_i;
            ins.read((char *)bin, sizeof(uint32_t) * 2 /* bin_id and n_chunk */);
            // allocate chunks
            bin->chunk = new index_chunk_t[bin->n_chunk];
            ins.read((char *)bin->chunk, sizeof(index_chunk_t) * bin->n_chunk);
        }

        // read n_intv
        ins.read((char *)&(ref->n_intv), sizeof(uint32_t));
        // allocate ioffset
        ref->ioffset = new uint64_t[ref->n_intv];
        // read ioffsets
        ins.read((char *)ref->ioffset, sizeof(uint64_t) * ref->n_intv);
    }
    return index;
}

//! Create parallelization regions from index
//! \param index The index object to a bam file
//! \param filesize The total file size of the bam file the index was created on
//! \return A vector of [ioffset_begin, ioffset_end) pairs, which can be used to parallelize compute jobs
std::vector<region> index_to_regions(const index_t& index, uint64_t filesize);

#endif
