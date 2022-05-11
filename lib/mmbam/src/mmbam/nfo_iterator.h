#ifndef NFO_ITERATOR_H
#define NFO_ITERATOR_H

#include <unistd.h>
#include <iostream>

//! A generic iterator for memory buffers with known offset
/*! Data structure like alignment reads in a decompressed bam buffer follows a 
 *  format of singly linked list. The location of the next buffer is known based
 *  on data stored in the current buffer. This iterator will automatically
 *  calcuate the correct address of the next record based on a known type of
 *  "block-size" of the current record, plus an optional fixed offset that are
 *  not included in the block-size field.
 *
 *  Note that this type is used to build other types of libmmbam, and is not intended
 *  to be used directly by client code
 *
 *  \tparam T The type of the record at each valid starting address
 *  \tparam OffT The type of the "block-size" data field
 *  \tparam byte_offset The number of bytes from the beginning of the record where the "block-size" field is located
 *  \tparam fixed_offset Any additional bytes the "block-size" field is not accounting for
 */
template<typename T, typename OffT, uint32_t byte_offset, uint32_t fixed_offset=0>
struct nfo_iterator : std::iterator<std::forward_iterator_tag, T> {
    const uint8_t* ptr;

    //! byte pointer to the next record
    const uint8_t* next_ptr() {
        const uint8_t* offset_ptr = ptr + byte_offset;
        auto total_offset = *reinterpret_cast<const OffT *>(offset_ptr) + fixed_offset;
        return ptr + total_offset;
    }

    //! construct nfo_iterator from a pointer
    nfo_iterator(const T* ptr) : ptr(reinterpret_cast<const uint8_t*>(ptr)) {}

    //! construct nfo_iterator from a std::vector of bytes
    nfo_iterator(const std::vector<uint8_t>& byte_buffer) : ptr(byte_buffer.data()) {}

    //! construct nfo_iterator from a std::vector of bytes, at a specified offset
    nfo_iterator(const std::vector<uint8_t>& byte_buffer, size_t offset) : ptr(byte_buffer.data() + offset) {}

    nfo_iterator operator++() { auto orig = *this; ptr = next_ptr(); return orig; }
    nfo_iterator operator++(int) { ptr = next_ptr(); return *this; }
    const T& operator*() { return *reinterpret_cast<const T*>(ptr); }
    const T* operator->() { return reinterpret_cast<const T*>(ptr); }
    off_t operator-(const nfo_iterator& rhs) { return ptr - rhs.ptr; }
    bool operator==(const nfo_iterator& rhs) { return ptr == rhs.ptr; }
    bool operator!=(const nfo_iterator& rhs) { return !(*this == rhs); }

    // casting nfo_iterator to a pointer returns the type the iterator is created with
    operator T const* () const { return reinterpret_cast<const T*>(ptr); }

    // casting nfo_iterator to a pointer of a arbitary type returns the pointer of that type
    template<typename PtrT>
    operator PtrT const* () const { return reinterpret_cast<const PtrT*>(ptr); }
};

#endif
