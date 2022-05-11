#ifndef MMBAM_MPILEUP_H
#define MMBAM_MPILEUP_H

#include <vector>
#include <list>
#include <cstring>
#include <functional>
#include "bam.h"
#include "index.h"

using pileup_filter_func_t = bool (*)(bam_rec_t* read);

//! Read-wise information at a piled up location
struct pileup_info_t {
    size_t buffer_id; /**< which buffer the read was from */
    size_t offset;    /**< where the read is located in the corresponding buffer */
    uint32_t qpos;    /**< the query position of the read that overlaps the piled up position */
    bool is_deletion; /**< whether the read contains a deletion at the piled up position */
};

#define PILEUP_MAX_DEPTH 4000

//! Information about a piled up location, produced by the mpileup engine and offered to the visitor function as a parameter
struct mpileup_t {
    int32_t ref_id;  /**< the numeric reference contig id of the piled up location */
    int32_t pos;     /**< the 0-based genomic coordinate of the piled up location */
    std::vector<const std::vector<uint8_t>*> reads_buffer; /**< The buffers (1 per file) of all reads */
    pileup_info_t *info; /**< additional, read-wise information at the piled up location, length= # of buffers * MAX_PILEUP_DEPTH */
    size_t *depth;       /**< number of reads overlapping this position, length = # of buffers */

    //! returns the correct info object
    //! \param b the index of the BAM file (or input buffer)
    //! \param d the index of the read within the specified buffer
    //! \return The info object corresponding to the specified buffer and read
    pileup_info_t& get_info(size_t b, size_t d) {
        return *(info + b * PILEUP_MAX_DEPTH + d);
    }

    //! returns the correct info object (const)
    //! \param b the index of the BAM file (or input buffer)
    //! \param d the index of the read within the specified buffer
    //! \return The info object corresponding to the specified buffer and read
    const pileup_info_t& get_info(size_t b, size_t d) const {
        return *(info + b * PILEUP_MAX_DEPTH + d);
    }

    void allocate(const std::vector<std::vector<uint8_t>>& buffers) {
        auto n_buffers = buffers.size();
        reads_buffer.resize(n_buffers);
        for(size_t i=0; i<n_buffers; i++) reads_buffer[i] = &buffers[i];
        info = new pileup_info_t[n_buffers * PILEUP_MAX_DEPTH];
        depth = new size_t[n_buffers];
    }

    void set_location(int32_t ref_id, int32_t pos, std::vector<std::vector<uint8_t>>& buffers)
    {
        memset(depth, 0, sizeof(size_t) * buffers.size());
        this->ref_id = ref_id;
        this->pos = pos;
    }

    void release() {
        if(info != nullptr) delete [] info;
        if(depth != nullptr) delete [] depth;
    }

};


//! Type alias of multiple input mfile
using mfiles_t = std::vector<std::reference_wrapper<const mfile_t::ptr_t>>;

//! Type alias of multiple index
using indices_t = std::vector<std::reference_wrapper<const index_t>>;

//! Multiple input pileup engine
//
//! \param mfiles A list of input memory mapped BAM files
//! \param indices A list of indices of the input BAM files, order-matched
//! \param ref_id The numerical reference contig id of the pileup region
//! \param pos_start The genomic coordinate of the start of the pileup region
//! \param pos_end The genomic coordinate of the end of the pileup region
//! \param predicate_func A functor, invoked with one parameter of a bam_rec_t struct const reference, and returns true if such read should be considered by the pileup engine
//! \param visitor_func A functor, invoked at each piled up location with one parameter of a mpileup_t struct const reference, and returns true if pileup engine should continue
void mpileup(mfiles_t mfiles,
        indices_t indices,
        uint32_t ref_id, int32_t pos_start, int32_t pos_end,
        const std::function<bool(const bam_rec_t&)>& predicate_func,
        const std::function<bool(const mpileup_t&)>& visitor_func);


#endif
