// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 @file real.h
 SYNOPSIS: we define and use "real" to be able to easily change
 the floating point precision, depending on the application.
 REAL_EPSILON is a lower limit of the precision achieved.
 */

#ifndef REAL_H
#define REAL_H


#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <new>

/**
 It is possible to select double or single precision throughout Cytosim
 by editing this file.
 
 Cytosim is faster in single precision, but the iterative solver used
 in Meca::solve() (conjugate-gradient) may fail in adverse conditions.
 Much of the code was not optimized for double precision.
 
 It is safer, and STRONGLY ADVISED therefore, to use double precision!
*/
#define REAL_IS_DOUBLE 1


#if REAL_IS_DOUBLE
   /// real is an alias to double
   typedef double real;
   constexpr real REAL_EPSILON = 128 * DBL_EPSILON;
#else
   /// real is an alias to float
   typedef float real;
   constexpr real REAL_EPSILON = 128 * FLT_EPSILON;
#endif


/// return a number greater or equal to 's' that is a multiple of 4
inline size_t chunk_real(size_t s)
{
    // align to 4 doubles (of size 8 bytes), hence 32 bytes
    constexpr size_t chunk = 32 / sizeof(real);
    // return a multiple of chunk greater than 's'
    // this bit trickery works if chunk is a pure power of 2
    return ( s + chunk - 1 ) & ~( chunk - 1 );
}

/// check memory alignement of a pointer
inline void check_alignment(void * ptr)
{
    uintptr_t a = ((uintptr_t)ptr & 63);
    if ( a )
        fprintf(stderr, "missaligned pointer %p (%lu)\n", ptr, a);
}


/// allocate a new array to hold `size` real scalars
/** The returned pointer is aligned to a 32 byte boundary */
inline real* new_real(size_t size)
{
    void* ptr = nullptr;
    // we align to 4 doubles (of size 8 bytes), hence 32 bytes
    if ( posix_memalign(&ptr, 32, size*sizeof(real)) )
        throw std::bad_alloc();
    //printf("new_real(%lu)  %lu\n", size, ((uintptr_t)ptr&63));
    return (real*)ptr;
}


/// release an array of reals allocated by `new_real`
inline void free_real(void * ptr)
{
    free(ptr);
}


/// copy `size` real scalars from `src` to `dst`
inline void copy_real(size_t size, real const* src, real * dst)
{
#if ( 0 )
    memcpy(dst, src, size*sizeof(real));
#else
    //#pragma vector unaligned
    for ( size_t u = 0; u < size; ++u )
        dst[u] = src[u];
#endif
}


/// set `size` values of the array `ptr` to 0 (zero).
inline void zero_real(size_t size, real * ptr)
{
#if ( 1 )
    memset(ptr, 0, size*sizeof(real));
#else
    #pragma vector unaligned
    for ( size_t u = 0; u < size; ++u )
        ptr[u] = 0.0;
#endif
}


/// square of the argument: `x * x`
inline real square(const real x) { return x * x; }

/// cube of the argument: `x * x * x`
inline real cube(const real x) { return x * x * x; }


#endif
