// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Monday 5 June 2018 was a very nice day in Strasbourg

#ifndef SIMD_H
#define SIMD_H

#define CHECK_VECTOR_ALIGNMENT 0

#include <cstdio>


//---------------------------------- SSE ---------------------------------------

#ifdef __SSE3__

#include <pmmintrin.h>

/// Vector of 2 doubles
typedef __m128d vec2;

inline void print(vec2 v, char const* s)
{
    printf("vec2 %s ( %5.2f %5.2f )\n", s, v[1], v[0]);
}

inline void print(vec2 v, vec2 w, char const* s)
{
    printf("vec2 %s ( %5.2f %5.2f )( %5.2f %5.2f )\n", s, v[1], v[0], w[1], w[0]);
}

/*
inline void print8(__m128i v, char const* s)
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

inline vec2 load1(double const* a)           { return _mm_load_sd(a); }
#if CHECK_VECTOR_ALIGNMENT
inline vec2 load2(double const* a)
{
    if ((uintptr_t)a & 15)
        fprintf(stderr, "unaligned __m128d load2 %p\n", a);
    return _mm_load_pd(a);
}
#else
inline vec2 load2(double const* a)           { return _mm_load_pd(a); }
#endif

// unaligned load
inline vec2 loadu2(double const* a)          { return _mm_loadu_pd(a); }

// load 1 double and duplicate
inline vec2 loaddup2(double const* a)        { return _mm_load1_pd(a); }

inline vec2 loadhi2(vec2 a, double const* b) { return _mm_loadh_pd(a,b); }
inline vec2 loadlo2(vec2 a, double const* b) { return _mm_loadl_pd(a,b); }

inline void store2(double* a, vec2 b)        { _mm_store_pd(a, b); }
inline void store12(double* a, vec2 b)       { _mm_store1_pd(a, b); }
inline void storelo(double* a, vec2 b)       { _mm_store_sd(a, b); }
inline void storeu2(double* a, vec2 b)       { _mm_storeu_pd(a, b); }

inline vec2 movedup2(vec2 a)                 { return _mm_movedup_pd(a); }

inline vec2 mul2(vec2 a, vec2 b)             { return _mm_mul_pd(a,b); }
inline vec2 div2(vec2 a, vec2 b)             { return _mm_div_pd(a,b); }
inline vec2 add2(vec2 a, vec2 b)             { return _mm_add_pd(a,b); }
inline vec2 sub2(vec2 a, vec2 b)             { return _mm_sub_pd(a,b); }
inline vec2 hadd2(vec2 a, vec2 b)            { return _mm_hadd_pd(a,b); }

inline vec2 sqrt2(vec2 a)                    { return _mm_sqrt_pd(a); }
inline vec2 max2(vec2 a, vec2 b)             { return _mm_max_pd(a,b); }
inline vec2 min2(vec2 a, vec2 b)             { return _mm_min_pd(a,b); }
inline vec2 and2(vec2 a, vec2 b)             { return _mm_and_pd(a,b); }
inline vec2 andnot2(vec2 a, vec2 b)          { return _mm_andnot_pd(a,b); }
inline vec2 abs2(vec2 a)                     { return _mm_andnot_pd(_mm_set1_pd(-0.0), a); }

inline vec2 setr2(double a, double b)        { return _mm_setr_pd(a,b); }
inline vec2 set2(double a, double b)         { return _mm_set_pd(a, b); }
inline vec2 set2(double a)                   { return _mm_set1_pd(a); }
inline vec2 setzero2()                       { return _mm_setzero_pd(); }

inline vec2 unpacklo2(vec2 a, vec2 b)        { return _mm_unpacklo_pd(a,b); }
inline vec2 unpackhi2(vec2 a, vec2 b)        { return _mm_unpackhi_pd(a,b); }

#define shuffle2(a,b,c)   _mm_shuffle_pd(a,b,c)
#define blend2(a,b,c)     _mm_blend_pd(a,b,c)
#define blendv2(a,b,c)    _mm_blendv_pd(a,b,c)
#define cmp2(a,b,c)       _mm_cmp_pd(a,b,c)

/// returns the sum of the elements, broadcasted
inline vec2 esum(vec2 v)
{
    return add2(v, shuffle2(v, v, 0b01));
}

/// returns the dot product of two vectors, broadcasted
inline vec2 dot2(vec2 a, vec2 b)
{
    vec2 p = mul2(a, b);
    return add2(p, shuffle2(p, p, 0b01));
}

/// square of vector norm, broadcasted
inline vec2 normsqr2(vec2 vec)
{
    vec2 p = mul2(vec, vec);
    return add2(p, shuffle2(p, p, 0b01));
}

/// normalize vector
inline vec2 normalize2(vec2 vec)
{
    vec2 p = mul2(vec, vec);
    vec2 s = add2(p, shuffle2(p, p, 0b01));
    return div2(vec, sqrt2(s));
}

/// normalize vector to 'n'
inline vec2 normalize2(vec2 vec, double n)
{
    vec2 p = mul2(vec, vec);
    vec2 s = add2(p, shuffle2(p, p, 0b01));
    return mul2(vec, div2(set2(n), sqrt2(s)));
}

typedef __m128  vec4f;

inline void print(vec4f v, char const* s)
{
    printf("vec4f %s ( %5.2f %5.2f %5.2f %5.2f )\n", s, v[3], v[2], v[1], v[0]);
}

inline vec4f load4f(float const* a)      { return _mm_load_ps(a); }
inline void store4f(float* a, vec4f b)   { return _mm_store_ps(a, b); }
inline vec4f max4f(vec4f a, vec4f b)     { return _mm_max_ps(a,b); }
inline vec4f min4f(vec4f a, vec4f b)     { return _mm_min_ps(a,b); }
inline vec4f and4f(vec4f a, vec4f b)     { return _mm_and_ps(a,b); }
inline vec4f andnot4f(vec4f a, vec4f b)  { return _mm_andnot_ps(a,b); }
inline vec4f abs4f(vec4f a)              { return _mm_andnot_ps(_mm_set1_ps(-0.0), a); }
#define permute4f(a,b)       _mm_permute_ps(a,b)       // same as shuffle2(a,a,b)

#endif

//---------------------------------- AVX ---------------------------------------

#ifdef __AVX__

#include <immintrin.h>

/// Vector of 4 doubles
typedef __m256d vec4;

inline void print(vec4 v, char const* s)
{
    printf("vec4 %s ( %5.2f %5.2f %5.2f %5.2f )\n", s, v[3], v[2], v[1], v[0]);
}

inline void print(vec4 v, vec4 w, char const* s)
{
    printf("vec4 %s ( %5.2f %5.2f %5.2f %5.2f )( %5.2f %5.2f %5.2f %5.2f )\n",
           s, v[3], v[2], v[1], v[0], w[3], w[2], w[1], w[0]);
}

/*
inline void print4(__m128i v, char const* s)
{
    uint32_t a = _mm_extract_epi32(v, 0);
    uint32_t b = _mm_extract_epi32(v, 1);
    uint32_t c = _mm_extract_epi32(v, 2);
    uint32_t d = _mm_extract_epi32(v, 3);
    printf("veci %s ( %5i %5i %5i %5i )\n", s, d, c, b, a);
}
*/

#define set64x(a,b,c,d)     _mm256_setr_epi64x(a,b,c,d)

constexpr __m256i msk3 = {-1,-1,-1,0};

//#define load3(a)            blend4(cast4(load2(a)), _mm256_broadcast_sd(a+2), 0b0100)
//#define store3(a,b)         storeu2(a, getlo(b)); storelo(Y+ii+2, gethi(z));

//inline vec4  load3(double const* a)  { return _mm256_loadu_pd(a); }
inline vec4 load3(double const* a)     { return _mm256_maskload_pd(a, msk3); }

/// load 4 values, and zeros out the upper one
inline vec4 load4z(double const* a)    { return _mm256_blend_pd(_mm256_loadu_pd(a), _mm256_setzero_pd(), 0b1000); }

#if CHECK_VECTOR_ALIGNMENT
inline vec4 load4(double const* a)
{
    if ((uintptr_t)a & 31)
        fprintf(stderr, "unaligned __m256d load4 %p (%lu)\n", a, (uintptr_t)a & 31);
    return _mm256_load_pd(a);
}
#else
inline vec4 load4(double const* a)              { return _mm256_load_pd(a); }
#endif
inline vec4 loadu4(double const* a)             { return _mm256_loadu_pd(a); }


inline void store3(double* a, vec4 b)           { _mm256_maskstore_pd(a, msk3, b); }
inline void store4(double* a, vec4 b)           { _mm256_store_pd(a,b); }
inline void storeu4(double* a, vec4 b)          { _mm256_storeu_pd(a,b); }

inline vec4 maskload4(double const* a, __m256i b)     { return _mm256_maskload_pd(a,b); }
inline void maskstore4(double* a, __m256i b, vec4 c)  { _mm256_maskstore_pd(a,b,c); }

inline vec4 setr4(double a, double b, double c, double d) { return _mm256_setr_pd(a,b,c,d); }
inline vec4 set4(double a, double b, double c, double d)  { return _mm256_set_pd(a,b,c,d); }

inline vec4 set4(double a)               { return _mm256_set1_pd(a); }
inline vec4 setzero4()                   { return _mm256_setzero_pd(); }

inline vec4 duplo4(vec4 a)               { return _mm256_movedup_pd(a); }
inline vec4 duphi4(vec4 a)               { return _mm256_permute_pd(a,15); }

/// load one double into all 4 positions
inline vec4 broadcast1(double const* a)  { return _mm256_broadcast_sd(a); }
/// load two double and duplicate: X, Y, X, Y
inline vec4 broadcast2(double const* a)  { return _mm256_broadcast_pd((__m128d const*)a); }

inline vec2 getlo(vec4 a)                { return _mm256_castpd256_pd128(a); }
inline vec2 gethi(vec4 a)                { return _mm256_extractf128_pd(a,1); }
inline vec4 cast4(vec2 a)                { return _mm256_castpd128_pd256(a); }
inline vec2 cast2(vec4 a)                { return _mm256_castpd256_pd128(a); }

inline vec4 mul4(vec4 a, vec4 b)         { return _mm256_mul_pd(a,b); }
inline vec4 div4(vec4 a, vec4 b)         { return _mm256_div_pd(a,b); }
inline vec4 add4(vec4 a, vec4 b)         { return _mm256_add_pd(a,b); }
inline vec4 sub4(vec4 a, vec4 b)         { return _mm256_sub_pd(a,b); }
inline vec4 hadd4(vec4 a, vec4 b)        { return _mm256_hadd_pd(a,b); }

inline vec4 sqrt4(vec4 a)                { return _mm256_sqrt_pd(a); }
inline vec4 max4(vec4 a, vec4 b)         { return _mm256_max_pd(a,b); }
inline vec4 min4(vec4 a, vec4 b)         { return _mm256_min_pd(a,b); }
inline vec4 and4(vec4 a, vec4 b)         { return _mm256_and_pd(a,b); }
inline vec4 andnot4(vec4 a, vec4 b)      { return _mm256_andnot_pd(a,b); }
inline vec4 abs4(vec4 a)                 { return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a); }

inline vec4 unpacklo4(vec4 a, vec4 b)    { return _mm256_unpacklo_pd(a,b); }
inline vec4 unpackhi4(vec4 a, vec4 b)    { return _mm256_unpackhi_pd(a,b); }

/* Unused functions:
 inline vec4 loadu22(double const* a, double const* b) { return _mm256_loadu2_m128d(a,b); }
 inline void store22(double* a, double* b, vec4 c) { return _mm256_storeu2_m128d(a,b,c); }
 */

#define insertf128(a,b,c)   _mm256_insertf128_pd(a,b,c)
#define permute4(a,b)       _mm256_permute_pd(a,b)
#define permute2(a,b)       _mm_permute_pd(a,b)       // same as shuffle2(a,a,b)
#define permute2f128(a,b,c) _mm256_permute2f128_pd(a,b,c)
#define shuffle4(a,b,c)     _mm256_shuffle_pd(a,b,c)
#define blend4(a,b,mask)    _mm256_blend_pd(a,b,mask)
#define blendv4(a,b,mask)   _mm256_blendv_pd(a,b,mask)
#define cmp4(a,b,c)         _mm256_cmp_pd(a,b,c)


/// concatenate two vec2 into a vec4
inline vec4 cat4(vec2 h, vec2 l) { return _mm256_insertf128_pd(_mm256_castpd128_pd256(l), h, 1); }
inline vec4 cat4(vec2 h, vec4 l) { return _mm256_insertf128_pd(l, h, 1); }

//inline vec4 cat4(vec2 h, vec2 l) { return _mm256_set_m128d(h, l); }
//#define cat4(h, l)           _mm256_set_m128d(h, l)


/// returns the sum of the elements, broadcasted
inline vec4 esum(vec4 v)
{
    vec4 s = add4(v, permute2f128(v, v, 0x01));
    return add4(s, permute4(s, 0b0101));
}

/// returns the dot product of two vectors, broadcasted
inline vec4 dot4(vec4 a, vec4 b)
{
    vec4 m = mul4(a, b);
    vec4 s = add4(m, permute2f128(m, m, 0x01));
    return add4(s, permute4(s, 0b0101));
}

/// square of vector norm, broadcasted
inline vec4 normsqr4(vec4 vec)
{
    vec4 m = mul4(vec, vec);
    vec4 s = add4(m, permute2f128(m, m, 0x01));
    return add4(s, permute4(s, 0b0101));
}

/// normalize vector
inline vec4 normalize4(vec4 vec)
{
    vec4 m = mul4(vec, vec);
    vec4 s = add4(m, permute2f128(m, m, 0x01));
    m = add4(s, permute4(s, 0b0101));
    return div4(vec, sqrt4(m));
}

/// normalize vector to 'n'
inline vec4 normalize4(vec4 vec, double n)
{
    vec4 m = mul4(vec, vec);
    vec4 s = add4(m, permute2f128(m, m, 0x01));
    m = add4(s, permute4(s, 0b0101));
    return mul4(vec, div4(set4(n), sqrt4(m)));
}

#endif

//-------------------------- AVX Single Precision-------------------------------

#ifdef __AVX__

/// Vector of 8 floats
typedef __m256 vec8f;

inline void print(vec8f v, char const* x)
{
    printf("vec8f %s ( %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f )\n", x,
           v[7], v[6], v[5], v[4], v[3], v[2], v[1], v[0]);
}

inline vec8f load8f(float const* a)     { return _mm256_load_ps(a); }
inline void store8f(float* a, vec8f b)  { return _mm256_store_ps(a, b); }

inline vec8f setzero8f()                { return _mm256_setzero_ps(); }
inline vec8f set8f(float const& a)      { return _mm256_set1_ps(a); }
inline vec8f mulf(vec8f a, vec8f b)     { return _mm256_mul_ps(a,b); }
inline vec8f divf(vec8f a, vec8f b)     { return _mm256_div_ps(a,b); }
inline vec8f addf(vec8f a, vec8f b)     { return _mm256_add_ps(a,b); }
inline vec8f subf(vec8f a, vec8f b)     { return _mm256_sub_ps(a,b); }

inline vec8f max8f(vec8f a, vec8f b)    { return _mm256_max_ps(a,b); }
inline vec8f min8f(vec8f a, vec8f b)    { return _mm256_min_ps(a,b); }
inline vec8f and8f(vec8f a, vec8f b)    { return _mm256_and_ps(a,b); }
inline vec8f andnot8f(vec8f a, vec8f b) { return _mm256_andnot_ps(a,b); }
inline vec8f abs8f(vec8f a)             { return _mm256_andnot_ps(_mm256_set1_ps(-0.0), a); }

#define permute8f128(a,b,c)  _mm256_permute4f128_ps(a,b,c)

/// approximate inverse
inline vec8f rcpf(vec8f a)              { return _mm256_rcp_ps(a); }
/// approximate inverse square root
inline vec8f rsqrtf(vec8f a)            { return _mm256_rsqrt_ps(a); }

inline vec4  cvt4f(vec4f a)             { return _mm256_cvtps_pd(a); }
inline vec4f getlof(vec8f a)            { return _mm256_castps256_ps128(a); }
inline vec4f gethif(vec8f a)            { return _mm256_extractf128_ps(a,1); }

inline vec8f cvt8i(__m256i a)           { return _mm256_cvtepi32_ps(a); }

#define load8si(a)           _mm256_load_si256(a)
#define cmpf(a,b,c)          _mm256_cmp_ps(a,b,c)
#define permute2f128f(a,b,c) _mm256_permute2f128_ps(a,b,c)

#endif

//---------------------------------- AVX2 --------------------------------------

#ifdef __AVX2__

#define permute4x64(a,b)    _mm256_permute4x64_pd(a,b)
#define rotater4(a)         _mm256_castsi256_pd(_m256_alignr_epi8(a, a, 1));
#define rotatel4(a)         _mm256_castsi256_pd(_m256_alignr_epi8(a, a, 7));

// load a memory address = [ X Y ] into [ X X Y Y ]:
#define interleave4(a)      _mm256_permute4x64_pd(cast4(a), 0x50)

inline vec4 broadcast1(vec2 a)  { return _mm256_broadcastsd_pd(a); }


/// cross product of two 3D vectors ( X Y Z T )
inline vec4 cross4(vec4 a, vec4 b)
{
    vec4 a1 = permute4x64(a, 0xC9); // Y Z X T
    vec4 b1 = permute4x64(b, 0xC9); // Y Z X T
    return permute4x64(sub4(mul4(a,b1), mul4(a1,b)), 0xC9);
}

#else

#define interleave4(a)      permute4(permute2f128(cast4(a), cast4(a), 0x00), 0b1100)

#endif

//----------------------------------- FMA --------------------------------------

#ifdef __FMA__

inline vec2 fmadd1(vec2 a, vec2 b, vec2 c)  { return _mm_fmadd_sd(a,b,c); }  // a * b + c
inline vec2 fmsub1(vec2 a, vec2 b, vec2 c)  { return _mm_fmsub_sd(a,b,c); }  // a * b - c
inline vec2 fnmadd1(vec2 a, vec2 b, vec2 c) { return _mm_fnmadd_sd(a,b,c); } // c - a * b

inline vec2 fmadd2(vec2 a, vec2 b, vec2 c)  { return _mm_fmadd_pd(a,b,c); }
inline vec2 fmsub2(vec2 a, vec2 b, vec2 c)  { return _mm_fmsub_pd(a,b,c); }
inline vec2 fnmadd2(vec2 a, vec2 b, vec2 c) { return _mm_fnmadd_pd(a,b,c); }

inline vec4 fmadd4(vec4 a, vec4 b, vec4 c)  { return _mm256_fmadd_pd(a,b,c); }
inline vec4 fmsub4(vec4 a, vec4 b, vec4 c)  { return _mm256_fmsub_pd(a,b,c); }
inline vec4 fnmadd4(vec4 a, vec4 b, vec4 c) { return _mm256_fnmadd_pd(a,b,c); }

#else

// define erzatz functions
//#warning "Patching SIMD' Fused Multiply Add functions"

#ifdef __SSE3__
inline vec2 fmadd1(vec2 a, vec2 b, vec2 c)  { return _mm_add_sd(_mm_mul_sd(a,b), c); }
inline vec2 fmsub1(vec2 a, vec2 b, vec2 c)  { return _mm_sub_sd(_mm_mul_sd(a,b), c); }
inline vec2 fnmadd1(vec2 a, vec2 b, vec2 c) { return _mm_sub_sd(c, _mm_mul_sd(a,b)); }

inline vec2 fmadd2(vec2 a, vec2 b, vec2 c)  { return _mm_add_pd(_mm_mul_pd(a,b), c); }
inline vec2 fmsub2(vec2 a, vec2 b, vec2 c)  { return _mm_sub_pd(_mm_mul_pd(a,b), c); }
inline vec2 fnmadd2(vec2 a, vec2 b, vec2 c) { return _mm_sub_pd(c, _mm_mul_pd(a,b)); }
#endif

#ifdef __AVX__
inline vec4 fmadd4(vec4 a, vec4 b, vec4 c)  { return _mm256_add_pd(_mm256_mul_pd(a,b), c); }
inline vec4 fmsub4(vec4 a, vec4 b, vec4 c)  { return _mm256_sub_pd(_mm256_mul_pd(a,b), c); }
inline vec4 fnmadd4(vec4 a, vec4 b, vec4 c) { return _mm256_sub_pd(c, _mm256_mul_pd(a,b)); }
#endif

#endif

#endif

