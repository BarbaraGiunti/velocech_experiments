/*  Copyright 2013 IST Austria
    Contributed by: Hubert Wagner

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */


/////////////// include/phat/helpers/misc.h

#pragma once
// STL includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <queue>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <iterator>

// VS2008 and below unfortunately do not support stdint.h
#if defined(_MSC_VER)&& _MSC_VER < 1600
    typedef __int8 int8_t;
    typedef unsigned __int8 uint8_t;
    typedef __int16 int16_t;
    typedef unsigned __int16 uint16_t;
    typedef __int32 int32_t;
    typedef unsigned __int32 uint32_t;
    typedef __int64 int64_t;
    typedef unsigned __int64 uint64_t;
#else
    #include <stdint.h>
#endif

// basic types. index can be changed to int32_t to save memory on small instances
namespace phat {
    typedef int64_t index;
    typedef int8_t dimension;
    typedef std::vector< index > column;
}

// OpenMP (proxy) functions
#if defined _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
	inline void omp_set_num_threads( int ) {};
    #include <time.h>
    #define omp_get_wtime() (float)clock() / (float)CLOCKS_PER_SEC
#endif


/////////////// include/phat/helpers/thread_local_storage.h


// should ideally be equal to the cache line size of the CPU
#define PHAT_TLS_SPACING_FACTOR 64

// ThreadLocalStorage with some spacing to avoid "false sharing" (see wikipedia)
template< typename T > 
class thread_local_storage
{
public:

    thread_local_storage() : per_thread_storage( omp_get_max_threads() * PHAT_TLS_SPACING_FACTOR ) {};

    T& operator()() {
        return per_thread_storage[ omp_get_thread_num() * PHAT_TLS_SPACING_FACTOR ];
    }

    const T& operator()() const {
        return per_thread_storage[ omp_get_thread_num() * PHAT_TLS_SPACING_FACTOR ];
    }

    T& operator[]( int tid ) {
        return per_thread_storage[ tid * PHAT_TLS_SPACING_FACTOR ];
    }

    const T& operator[]( int tid ) const {
        return per_thread_storage[ tid * PHAT_TLS_SPACING_FACTOR ];
    }

protected:
    std::vector< T > per_thread_storage; 
};




/////////////// include/phat/representations/bit_tree_pivot_column.h





namespace phat {

    // This is a bitset indexed with a 64-ary tree. Each node in the index
    // has 64 bits; i-th bit says that the i-th subtree is non-empty.
    // Supports practically O(1), inplace, zero-allocation: insert, remove, max_element
    // and clear in O(number of ones in the bitset).
    // 'add_index' is still the real bottleneck in practice.
    class bit_tree_column
    {
    protected:

        size_t offset; // data[i + offset] = ith block of the data-bitset
        typedef uint64_t block_type;
        std::vector< block_type > data;

        
        size_t debruijn_magic_table[ 64 ];

        enum { block_size_in_bits = 64 };
        enum { block_shift = 6 };

        // Some magic: http://graphics.stanford.edu/~seander/bithacks.html
        // Gets the position of the rightmost bit of 'x'. 0 means the most significant bit. 
        // (-x)&x isolates the rightmost bit.
        // The whole method is much faster than calling log2i, and very comparable to using ScanBitForward/Reverse intrinsic,
        // which should be one CPU instruction, but is not portable.
        size_t rightmost_pos( const block_type value ) const {                
            return 64 - 1 - debruijn_magic_table[ ( (value & (-(int64_t)value) ) * 0x07EDD5E59A4E28C2 ) >> 58 ];
        }

    public:        

        void init( index num_cols ) {
            int64_t n = 1; // in case of overflow
            int64_t bottom_blocks_needed = ( num_cols + block_size_in_bits - 1 ) / block_size_in_bits;
            int64_t upper_blocks = 1;        

            // How many blocks/nodes of index needed to index the whole bitset?
            while( n * block_size_in_bits < bottom_blocks_needed ) {
                n *= block_size_in_bits;
                upper_blocks += n;
            }

            offset = upper_blocks;
            data.resize( upper_blocks + bottom_blocks_needed, 0 );

            std::size_t temp_array[ 64 ] = {
                    63,  0, 58,  1, 59, 47, 53,  2,
                    60, 39, 48, 27, 54, 33, 42,  3,
                    61, 51, 37, 40, 49, 18, 28, 20,
                    55, 30, 34, 11, 43, 14, 22,  4,
                    62, 57, 46, 52, 38, 26, 32, 41,
                    50, 36, 17, 19, 29, 10, 13, 21,
                    56, 45, 25, 31, 35, 16,  9, 12,
                    44, 24, 15,  8, 23,  7,  6,  5 };

            std::copy( &temp_array[ 0 ], &temp_array[ 64 ], &debruijn_magic_table[ 0 ] );
        }

        index get_max_index() const {
            if( !data[ 0 ] )
                return -1;

            size_t n = 0;
            size_t newn = 0;
            size_t index = 0;
            while( newn < data.size() ) {
                n = newn;
                index = rightmost_pos( data[ n ] );
                newn = ( n << block_shift ) + index + 1;
            } 

            return ( ( n - offset ) << block_shift ) + index;
        }

        bool is_empty() const {
            return data[ 0 ] == 0;
        }

        void add_index( const size_t entry ) {
            const block_type ONE = 1;
            const block_type block_modulo_mask = ( ONE << block_shift ) - 1;
            size_t index_in_level = entry >> block_shift;
            size_t address = index_in_level + offset;
            size_t index_in_block = entry & block_modulo_mask;

            block_type mask = ( ONE << ( block_size_in_bits - index_in_block - 1 ) );

            data[ address ] ^= mask;

            // Check if we reached the root. Also, if anyone else was in this block, we don't need to update the path up.
            while( address && !( data[ address ] & ~mask ) ) {
                index_in_block = index_in_level & block_modulo_mask;
                index_in_level >>= block_shift;
                --address;
                address >>= block_shift;
                mask = ( ONE << ( block_size_in_bits - index_in_block - 1 ) );
                data[ address ] ^= mask;
            }
        }

        void get_col_and_clear( column &out ) {
	    out.clear();
            index mx = this->get_max_index();
            while( mx != -1 ) {
                out.push_back( mx );
                add_index( mx );
                mx = this->get_max_index();
            }

            std::reverse( out.begin(), out.end() );
        }

	
	template<typename InputIterator>
	  void add_col(InputIterator begin, InputIterator end) {
	  for(InputIterator it=begin;it!=end;it++) {
	    add_index(*it);
	  }
	}
	

        void add_col(const column &col) {
            for( size_t i = 0; i < col.size(); ++i ){
                //std::cout << i << ",,, ";
                //std::cout << col.size() << std::endl;
                //std::cout << col[i] <<  ",, ";
                add_index(col[i]);
            }
        }

	void clear() {
            index mx = this->get_max_index();
            while( mx != -1 ) {
                add_index( mx );
                mx = this->get_max_index();
            }
        }

	index size() {
	    // INEFFICIENT!!
	    column col;
	    get_col(col);
	    return col.size();
	}
	
	void remove_max() {
            add_index( get_max_index() );
        }

        void set_col( const column& col ) {
            clear();
            add_col( col );
        }

        void get_col( column& col ) {
            get_col_and_clear( col );
            add_col( col );
        }
    };
}