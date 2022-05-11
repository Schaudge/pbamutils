#include "mmbam/mbgzf.h"
#include <vector>
#include <exception>
#include <iostream>
#include <functional>
#include <libdeflate.h>

#ifdef USE_ZLIB
#include <zlib.h>
static size_t bgzf_decompress(const uint8_t *src, size_t src_len, uint8_t *dest, size_t dest_len);
#endif

static size_t libdeflate_decompress(const uint8_t *src, size_t src_len, uint8_t *dest, size_t dest_len);

std::vector<uint8_t> bgzf_inflate(const bgzf_block_t& block) {
    // check isize
    uint32_t isize = bgzf_isize(block);
    if(isize == 0) return std::vector<uint8_t>();

    std::vector<uint8_t> buffer(isize);
    //auto bytes_inf = bgzf_decompress((const uint8_t *)&block, block.bsize + 1, &buffer[0], isize);
    auto bytes_inf = libdeflate_decompress((const uint8_t *)&block, block.bsize + 1, &buffer[0], isize);

    if(isize != bytes_inf) throw std::runtime_error("isize mismatch");

#ifdef USE_ZLIB
    // check crc32
    uint32_t block_crc32 = *(uint32_t *)((uint8_t*)(&block) + (block.bsize +1) - 8);
    uint32_t compute_crc32 = crc32(0L, Z_NULL, 0);
    compute_crc32 = crc32(compute_crc32, &buffer[0], bytes_inf);
    if(compute_crc32 != block_crc32) throw std::runtime_error("crc32 mismatch");
#endif

    return buffer;
}

std::vector<uint8_t> bgzf_inflate_range(const uint8_t *src, const size_t src_len) {
    if(src_len == 0) return std::vector<uint8_t>();

    bgzf_iterator_t bgzf_head(reinterpret_cast<const bgzf_block_t*>(src));
    bgzf_iterator_t bgzf_it(reinterpret_cast<const bgzf_block_t*>(src));
    bgzf_iterator_t bgzf_last(reinterpret_cast<const bgzf_block_t*>(src + src_len));

    // calculate total isize
    size_t isize = 0;
    while(bgzf_it < bgzf_last) {
        auto this_isize = bgzf_isize(*bgzf_it);
        isize += this_isize;
        bgzf_it++;
    }


    std::vector<uint8_t> buffer(isize);

    // decompress in batches
    bgzf_it = bgzf_iterator_t(reinterpret_cast<const bgzf_block_t*>(src));
    size_t batch_isize = 0;
    unsigned int maxBufSize = - 1;

    size_t buffer_start = 0;
    while(bgzf_it < bgzf_last) {
        auto this_isize = bgzf_isize(*bgzf_it);
        batch_isize += this_isize;
        if(batch_isize > maxBufSize) {
            // decompress [bgzf_head, bgzf_it)
            //auto bytes_inf = bgzf_decompress(bgzf_head.ptr, bgzf_it.ptr - bgzf_head.ptr, &buffer[buffer_start], batch_isize - this_isize);
            auto bytes_inf = libdeflate_decompress(bgzf_head.ptr, bgzf_it.ptr - bgzf_head.ptr, &buffer[buffer_start], batch_isize - this_isize);
            // buffer is filled from buffer_start - buffer_start + isize;
            // next batch should fill starting from buffer_start + isize;
            if(bytes_inf != batch_isize - this_isize) throw std::runtime_error("batch isize mismatch");
            buffer_start = buffer_start + bytes_inf;
            bgzf_head = bgzf_it;
            batch_isize = this_isize;
        }
        bgzf_it++;
    }
    // one more fill
    //auto bytes_inf = bgzf_decompress(bgzf_head.ptr, bgzf_it.ptr - bgzf_head.ptr, &buffer[buffer_start], batch_isize);
    auto bytes_inf = libdeflate_decompress(bgzf_head.ptr, bgzf_it.ptr - bgzf_head.ptr, &buffer[buffer_start], batch_isize);
    if(bytes_inf != batch_isize) throw std::runtime_error("batch isize mismatch");
    
    if(buffer.size() != isize) throw std::runtime_error("isize mismatch");
    return buffer;
}

std::vector<uint8_t> bgzf_inflate_range_p(const uint8_t *src, const size_t src_len,
        const std::function<std::vector<uint8_t>(
            const uint8_t*,
            const size_t,
            const size_t,
            const std::vector<size_t>&,
            const std::vector<size_t>&,
            size_t (*inflate)(const uint8_t*, size_t, uint8_t*, size_t))>& parallel_inflate) {
    if(src_len == 0) return std::vector<uint8_t>();
    // iterate over bgzf blocks to identify
    // 1. total decompression size
    // 2. each block's starting address and length
    // 3. each block's decompression address and length
    bgzf_iterator_t bgzf_head(reinterpret_cast<const bgzf_block_t*>(src));
    bgzf_iterator_t bgzf_it(reinterpret_cast<const bgzf_block_t*>(src));
    bgzf_iterator_t bgzf_last(reinterpret_cast<const bgzf_block_t*>(src + src_len));

    std::vector<size_t> src_off_vector;
    std::vector<size_t> dest_off_vector;

    size_t dest_off = 0;


    while(bgzf_it < bgzf_last) {
        auto isize = bgzf_isize(*bgzf_it);
        src_off_vector.push_back(bgzf_it.ptr - src);
        dest_off_vector.push_back(dest_off);
        dest_off += isize;
        bgzf_it++;
    }

    auto dest_len = dest_off;

    //return parallel_inflate(src, src_len, dest_len, src_off_vector, dest_off_vector, bgzf_decompress);
    return parallel_inflate(src, src_len, dest_len, src_off_vector, dest_off_vector, libdeflate_decompress);
}

#ifdef USE_ZLIB
static size_t bgzf_decompress(const uint8_t *src, size_t src_len, uint8_t *dest, size_t dest_len) {
    z_stream strm = {0};
    int ret;

    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    strm.opaque = Z_NULL;

    strm.next_in = const_cast<Bytef*>(src);
    strm.avail_in = src_len;
    strm.next_out = (Bytef*)(dest);
    strm.avail_out = dest_len;

    if(inflateInit2(&strm, 16 + MAX_WBITS) != Z_OK) throw std::runtime_error("inflateInit2 failed");

    // decompress continuous gzip blocks
    while(strm.avail_in > 0) {
        ret = inflate(&strm, Z_FINISH);

        if(ret == Z_STREAM_END && strm.avail_in > 0) {
            if(inflateEnd(&strm) != Z_OK) throw std::runtime_error("inflateEnd failed");
            if(inflateInit2(&strm, 16 + MAX_WBITS) != Z_OK) throw std::runtime_error("inflateInit2 failed");
        }
        if(ret != Z_STREAM_END) {
            std::cerr<<"inflate failed, ret="<<ret<<std::endl;
            std::cerr<<"avail_in="<<strm.avail_in<<", avail_out="<<strm.avail_out<<std::endl;
            throw std::runtime_error("inflate failed");
        }
    }

    if(inflateEnd(&strm) != Z_OK) throw std::runtime_error("inflateEnd failed");

    return dest_len - strm.avail_out;
}
#endif

static size_t libdeflate_decompress(const uint8_t *src, size_t src_len, uint8_t *dest, size_t dest_len) {
    size_t actual_in_nbytes;
    size_t actual_out_nbytes;

    struct libdeflate_decompressor *d;
    d = libdeflate_alloc_decompressor();

    auto ret = libdeflate_gzip_decompress_ex(d, src, src_len, dest, dest_len, &actual_in_nbytes, &actual_out_nbytes);
    libdeflate_free_decompressor(d);

    return actual_out_nbytes;
}
