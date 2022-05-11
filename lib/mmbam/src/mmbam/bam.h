#ifndef BAM_H
#define BAM_H

/****************************************************************
 * BAM file structure APIs implemented on memory mapped file
*****************************************************************/

#include <stdint.h>
#include <vector>
#include "nfo_iterator.h"
#include "mfile.h"
#include "index.h"

#define BYTEREF(b) ((const uint8_t*)(b))
#define BAMREF(b)  ((const bam_rec_t *)(b))
#define BAM_NEXT(b) (BAMREF(BYTEREF(b) + (b)->block_size + 4))
#define READ_IN_REGION(b, rbegin, rend) ((b)->pos < (rend) && (b)->pos + bam_query_length(b) > (rbegin))

#define READ_PAIRED 0x1
#define READ_MAPPED_PROPER_PAIR 0x2
#define READ_UNMAPPED 0x4
#define MATE_UNMAPPED 0x8
#define READ_REVERSE_STRAND 0x10
#define MATE_REVERSE_STRAND 0x20
#define READ_FIRST_IN_PAIR 0x40
#define READ_SECOND_IN_PAIR 0x80
#define READ_SECONDARY_ALIGNMENT 0x100
#define READ_FAILED_QUALITY_CHECKS 0x200
#define READ_PCR_OPTICAL_DUPLICATE 0x400
#define READ_SUPPLEMENTARY 0x800
#define READ_PRIMARY 0x900

//! This struct mirrors the byte layout of uncompressed bam header
struct bam_header_t {
    char magic[4];      // BAM\1
    uint32_t l_text;    // Length of header text
    uint8_t vardata[];  // Data of variable length

    // char text[l_text];   // Header text
    // uint32_t n_ref;      // Number of reference sequences
    // struct[] {
    //   uint32_t l_name;   // Length of the reference name
    //   char name[l_name]; // Reference name
    //   uint32_t l_ref;    // Length of the reference sequence
    // }
};

//! this struct mirrors the byte layout of uncompressed bam records
struct bam_rec_t {
    uint32_t block_size;    /**< Total length of alignment record except this field */
    int32_t  ref_id;        /**< Reference sequence id; -1 if unmapped */
    int32_t  pos;           /**< 0-based left-most position */
    uint8_t  l_read_name;   /**< Length of read_name */
    uint8_t  mapq;          /**< Mapping quality */
    uint16_t bin;           /**< BAI index bin */
    uint16_t n_cigar_op;    /**< Number of operations in CIGAR */
    uint16_t flag;          /**< Bitwise flags */
    uint32_t l_seq;         /**< Length of SEQ */
    int32_t  next_ref_id;   /**< Reference id of the next segment */
    int32_t  next_pos;      /**< 0-based left-most position of the next segment */
    int32_t  tlen;          /**< Template length */
    uint8_t  vardata[];     /**< Data of variable length */
    
    //  char[l_read_name]       // Read name
    //  uint32_t[n_cigar_op]    // CIGAR
    //  uint8_t[(l_seq + 1)/2]  // 4-bit encoded read seq
    //  char[l_seq]             // Phred-scale base quality

    //  struct[] {              // TAGS
    //      char[2] tag;        // Two-character tag
    //      char val_type;      // Value type
    //      value;              // Value (type by val_type)
};

using bam_iterator = nfo_iterator<bam_rec_t, uint32_t, 0, sizeof(uint32_t)>;

//! check if the given buffer contains a complete bam header record
bool bam_buffer_contains_header(const std::vector<uint8_t>& buffer); 

//! loads a block of bam records from a region
//! \param mfile The memory mapped bam file
//! \param ioffst_first the ioffset of the beginning of the region. ioffset is the virtual offset defined per SAM file specification, section 4.1.1
//! \param ioffset_last the ioffset of the end of the region.
//! \return a byte vector of all bam files contained in the region
std::vector<uint8_t> bam_load_block(const mfile_t::ptr_t& mfile, uint64_t ioffset_first, uint64_t ioffset_last);

//! returns the number of complete bam records found in the given buffer
size_t bam_count_records(const std::vector<uint8_t>& buffer);

//! returns the total length of the query string of the bam record
int16_t bam_query_length(const bam_rec_t* b);

//! returns the read name of the bam record
inline std::string bam_read_name(const bam_rec_t* b) {
    return std::string(b->vardata, b->vardata + b->l_read_name - 1);
}

//! returns the pointer to the cigar operations array of the bam record
inline const uint32_t* bam_cigar_ptr(const bam_rec_t* b) {
    return reinterpret_cast<const uint32_t *>(b->vardata
        + sizeof(char) * b->l_read_name);
}

//! returns the pointer to the query sequence string of the bam record
inline const uint8_t* bam_seq_ptr(const bam_rec_t* b) {
    return b->vardata
        + sizeof(char) * b->l_read_name
        + sizeof(uint32_t) * b->n_cigar_op;
}

//! returns the pointer to the base quality array of the bam record
inline const char* bam_bqual_ptr(const bam_rec_t* b) {
    return reinterpret_cast<const char *>(b->vardata
            + sizeof(char) * b->l_read_name
            + sizeof(uint32_t) * b->n_cigar_op
            + sizeof(uint8_t) * ((b->l_seq + 1) / 2));
}

//! extract the sequence from a 4-bit compact representation sequence array at given position
inline uint8_t bam_unpack_base(const uint8_t* seq, uint32_t qpos) {
    if(qpos % 2 == 0) // even base, higher 4 bits
        return seq[qpos / 2] >> 4;
    else // odd base, lower 4 bits
        return seq[qpos / 2] & 0xF;
}

const uint8_t bam_a_lo = 0x01;
const uint8_t bam_c_lo = 0x02;
const uint8_t bam_g_lo = 0x04;
const uint8_t bam_t_lo = 0x08;
const uint8_t bam_a_hi = 0x10;
const uint8_t bam_c_hi = 0x20;
const uint8_t bam_g_hi = 0x40;
const uint8_t bam_t_hi = 0x80;



inline char byte2base_lo(const uint8_t base) {
    return base & bam_a_lo ? 'A' :
        base & bam_c_lo ? 'C' :
        base & bam_g_lo ? 'G' :
        base & bam_t_lo ? 'T' : 'N';
};

inline char byte2base_hi(const uint8_t base){
    return base & bam_a_hi ? 'A' :
        base & bam_c_hi ? 'C' :
        base & bam_g_hi ? 'G' :
        base & bam_t_hi ? 'T' : 'N';
};

//! Helper function to load all bam records from a specific genomic region
//
//! \param mfile The memory mapped BAM file
//! \param index The index created on the BAM file
//! \param reg_id The numerical id of the reference contig
//! \param region_start The genomic coordinate of the start of the region
//! \param region_end The genomic coordinate of the end of the region
//! \return A byte vector containing the bam records within the specified region
std::vector<uint8_t>
bam_load_region(const mfile_t::ptr_t& mfile,
                const index_t& index,
                int32_t ref_id, int32_t region_start, int32_t region_end);

#endif
