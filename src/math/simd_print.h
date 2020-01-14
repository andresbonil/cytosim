// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SIMD_PRINT_H
#define SIMD_PRINT_H

// Functions to print SIMD vectors

#include <cstdio>

//---------------------------------- SSE ---------------------------------------

#ifdef __SSE3__

#include <pmmintrin.h>

/// print SIMD vector of 2 doubles
inline void show(vec2 v, char const* s)
{
    printf("vec2 %s ( %5.2f %5.2f )\n", s, v[1], v[0]);
}

/// print two SIMD vector
inline void show(vec2 v, vec2 w, char const* s)
{
    printf("vec2 %s ( %5.2f %5.2f )( %5.2f %5.2f )\n", s, v[1], v[0], w[1], w[0]);
}

/*
inline void show8(__m128i v, char const* s)
{
    uint16_t a = _mm_extract_epi16(v, 0);
    uint16_t b = _mm_extract_epi16(v, 1);
    uint16_t c = _mm_extract_epi16(v, 2);
    uint16_t d = _mm_extract_epi16(v, 3);
    uint16_t e = _mm_extract_epi16(v, 5);
    uint16_t f = _mm_extract_epi16(v, 6);
    uint16_t g = _mm_extract_epi16(v, 7);
    uint16_t h = _mm_extract_epi16(v, 8);
    printf("veci %s ( %3i %3i %3i %3i %3i %3i %3i %3i )\n", s, h, g, f, e, d, c, b, a);
}
*/

/// print SIMD vector of 4 floats
inline void show(vec4f v, char const* s)
{
    printf("vec4f %s ( %5.2f %5.2f %5.2f %5.2f )\n", s, v[3], v[2], v[1], v[0]);
}

#endif

//---------------------------------- AVX ---------------------------------------

#ifdef __AVX__

/// print SIMD vector of 4 doubles
inline void show(vec4 v, char const* s)
{
    printf("vec4 %s ( %5.2f %5.2f %5.2f %5.2f )\n", s, v[3], v[2], v[1], v[0]);
}

/// print two SIMD vector of 4 doubles
inline void show(vec4 v, vec4 w, char const* s)
{
    printf("vec4 %s ( %5.2f %5.2f %5.2f %5.2f )( %5.2f %5.2f %5.2f %5.2f )\n",
           s, v[3], v[2], v[1], v[0], w[3], w[2], w[1], w[0]);
}

/*
inline void show4(__m128i v, char const* s)
{
    uint32_t a = _mm_extract_epi32(v, 0);
    uint32_t b = _mm_extract_epi32(v, 1);
    uint32_t c = _mm_extract_epi32(v, 2);
    uint32_t d = _mm_extract_epi32(v, 3);
    printf("veci %s ( %5i %5i %5i %5i )\n", s, d, c, b, a);
}
*/

/// print SIMD vector of 8 floats
inline void show(vec8f v, char const* x)
{
    printf("vec8f %s ( %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f )\n", x,
           v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0]);
}

#endif

#endif
