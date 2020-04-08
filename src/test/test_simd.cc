// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/*
 Tests Intel's Streaming SIMD
 F. Nedelec, HD, July 2013
 
 To compile: c++ -O4 tictoc.cc -mavx test.cc
 To generate assembly: c++ -S test.cc
 */

#include <cstdio>
#include "tictoc.h"
#include <stdint.h>
#include "vecprint.h"
#include "simd.h"
#include "simd_print.h"


typedef double real;

const unsigned SIZE = 1<<14;
real x[SIZE], y[SIZE];

void init()
{
    for ( unsigned ii=0; ii<SIZE; ++ii )
    {
        x[ii] = 1.0/real(SIZE-ii);
        y[ii] = real(SIZE-ii);
    }
}

void show(size_t len, const real* vec)
{
    for ( size_t n = 0; n < len; ++n )
    {
        if ( n % 4 == 3 )
            printf("%+5.2f  ", vec[n]);
        else
            printf("%+5.2f ", vec[n]);
    }
    printf("\n");
}

real scalar()
{
    real d = 0;
    for ( unsigned ii=0; ii<SIZE; ++ii )
        d += x[ii] * y[ii];
    return d;
}

real vector2()
{
    vec2 s = setzero2();
    for ( unsigned ii=0; ii<SIZE; ii+=2 )
        s = add2(s, mul2( load2(x+ii), load2(y+ii) ));
    _mm_empty();
    
    return esum2(s)[0];
}

#ifdef __AVX__

real vector4()
{
    vec4 s = setzero4();
    for ( unsigned ii=0; ii<SIZE; ii+=4 )
        s = add4(s, mul4( load4(x+ii), load4(y+ii) ));
    _mm_empty();
    
    return esum4(s)[0];
}


real vectorU()
{
    vec4 v0 = setzero4();
    vec4 v1 = setzero4();
    vec4 v2 = setzero4();
    vec4 v3 = setzero4();
    
    for ( unsigned ii=0; ii<SIZE; ii+=16 )
    {
        v0 = add4(v0, mul4( load4(x+ii   ), load4(y+ii   ) ));
        v1 = add4(v1, mul4( load4(x+ii+4 ), load4(y+ii+4 ) ));
        v2 = add4(v2, mul4( load4(x+ii+8 ), load4(y+ii+8 ) ));
        v3 = add4(v3, mul4( load4(x+ii+12), load4(y+ii+12) ));
    }
    
    vec4 s = add4(add4(v0, v1), add4(v2, v3));
    _mm_empty();
    
    return esum4(s)[0];
}

#endif

void run(real (*func)(), const char name[])
{
    const int rep = 1<<14;
    real a = 0, b = 0, c = 0, d = 0;
    real e = 0, f = 0, g = 0, h = 0;
    fprintf(stderr, "%s:  ", name);
    TicToc::tic();
    for ( int ii=0; ii<rep; ++ii )
    {
        a = (*func)();
        b = (*func)();
        c = (*func)();
        d = (*func)();
        e = (*func)();
        f = (*func)();
        g = (*func)();
        h = (*func)();
    }
    
    real s = ((a + b) + (c + d)) + ((e + f) + (g + h));
    double ms = TicToc::toc();
    fprintf(stderr, " %f :  %.0f ms\n", s, ms);
}


void test_swapSSE()
{
    vec2 a = set2(0, 1);
    vec2 b = set2(2, 3);
    show(a, "a");
    show(b, "b");
    
#ifdef __AVX__
    show(permute2(b,0b00), "permute 0b00");
    show(permute2(b,0b01), "permute 0b01");
    show(permute2(b,0b10), "permute 0b10");
    show(permute2(b,0b11), "permute 0b11");
#endif
    
    show(shuffle2(a,b,0b00), "0b00");
    show(shuffle2(a,b,0b01), "0b01");
    show(shuffle2(a,b,0b10), "0b10");
    show(shuffle2(a,b,0b11), "0b11");
    show(unpacklo2(a,b), "unpacklo");
    show(unpackhi2(a,b), "unpackhi");
    show(movedup2(a),  "movedup(a)");
    show(movedup2(b),  "movedup(b)");
}

#ifdef __AVX__

/**
 This will convert src = { XYZ XYZ XYZ XYZ }
 into dst = { XYZ? XYZ? XYZ? XYZ? }
 */
inline void deswizzle4(real const* src, real * dst)
{
    vec4 s0 = load4(src);
    vec4 s1 = load4(src+4);
    vec4 s2 = load4(src+8);
    store4(dst, s0);
    
    vec4 d1 = permute2f128(s0, s1, 0x21);
    d1 = shuffle4(d1, s1, 0b0101);
    store4(dst+4 , d1);
    
    vec4 d2 = blend4(s1, s2, 0b0011);
    d2 = permute2f128(d2, d2, 0b0001);
    store4(dst+8 , d2);
    
    vec4 d3 = permute2f128(s2, s2, 0x01);
    d3 = shuffle4(s2, d3, 0b0101);
    store4(dst+12, d3);
}


void test_deswizzle()
{
    const int SIZE = 1024;
    real a[SIZE*12];
    real b[SIZE*16];
    
    for ( int n = 0; n < 12*SIZE; ++n )
        a[n] = 1 + n % 12;
    
#pragma ivdep
    for ( int n = 0; n < SIZE; ++n )
        deswizzle4(a+12*n, b+16*n);
    
    for ( int n = 0; n < 4*SIZE; ++n )
        b[4*n+3] = 0;
    
    show(12, a);
    show(16, b);
}

void test_cat()
{
    printf("------ test_cat\n");
    vec2 y{1.0, 2.0};
    vec2 x{3.0, 4.0};
    
    show(cat4(x, y), "cat4(x, y)");
    
    x = set2(1.0, 2.0);
    y = set2(3.0, 4.0);
    
    show(cat4(x, y), "cat4(x, y)");
}

void test_load()
{
    printf("------ test_load\n");
    double mem[4] = { 1, 2, 3, 4 };
    vec4 x{1.0, 2.0, 3.0, 4.0};
    show(x, "set ");
    
    vec4 y = load4(mem);
    show(y, "load");
    
    show(cast4(load1(mem)), "cast4(load1)");
    show(cast4(load2(mem)), "cast4(load2)");

    show(load3(mem), "load3");
    show(cat4(load1(mem+2), load2(mem)), "cat4(load1, load2)");

    vec4 t = blend4(load4(mem), setzero4(), 0b1000);
    show(t, "blend(load4, zero)");


    vec2 u = load2(mem);
    vec4 a = broadcast2(mem);
    vec4 n = permute4(a, 0b1100);

    show(u, "src");
    show(n, "permute4(broadcast2)");
    show(interleave4(u), "interleave4");
    show(cat4(u,u), "cat4");
    //show(_mm256_set_m128d(u,u),    "set_m128d");
    show(insertf128(cast4(u),u,1), "insertf128");
    show(permute4(cat4(u,u), 0b1100), "permute4(cat4, 0b1100)");
    
    vec4 p{-1.0, -2.0, -3.0, -4.0};
    _mm_storeu_pd(mem, _mm256_castpd256_pd128(p));
    _mm_store_sd(mem+2,_mm256_extractf128_pd(p,1));
    show(load4(mem), "load(store3)");
}


void test_broadcast()
{
    printf("------ test_broadcast\n");
    double mem[4] = { 1, 2, 3, 4 };
    
    // using 4 loads
    show(broadcast1(mem  ), "mem[0]");
    show(broadcast1(mem+1), "mem[1]");
    show(broadcast1(mem+2), "mem[2]");
    show(broadcast1(mem+3), "mem[3]");
    
    // using 2 loads
    vec4 xyxy = broadcast2(mem);
    vec4 ztzt = broadcast2(mem+2);
    show(duplo4(xyxy), "x");
    show(duphi4(xyxy), "y");
    show(duplo4(ztzt), "z");
    show(duphi4(ztzt), "t");
    
    // using 1 load
    vec4 xyzt = load4(mem);
    xyxy = permute2f128(xyzt, xyzt, 0x00);
    ztzt = permute2f128(xyzt, xyzt, 0x11);
    show(duplo4(xyxy), "X");
    show(duphi4(xyxy), "T");
    show(duplo4(ztzt), "Z");
    show(duphi4(ztzt), "T");
}

__m256i make_mask2(long i)
{
    switch( i )
    {
        case  0: return _mm256_setr_epi64x( 0, 0, 0, 0);
        case  1: return _mm256_setr_epi64x(-1, 0, 0, 0);
        case  2: return _mm256_setr_epi64x(-1,-1, 0, 0);
        case  3: return _mm256_setr_epi64x(-1,-1,-1, 0);
        default: return _mm256_setr_epi64x(-1,-1,-1,-1);
    }
}

__m256i make_mask(long i)
{
    vec4 v{0.5, 1.5, 2.5, 3.5};
    return _mm256_castpd_si256(cmp4(v, set4(i), _CMP_LT_OQ));
}

void test_store()
{
    printf("------ test_store\n");
    double mem[4] = { 0, 0, 0, 0 };
    vec4 x{1.0, 2.0, 3.0, 4.0};
    show(x, "value");
    for ( int i = 0; i < 5; ++i )
    {
        __m256i k = make_mask(i);
        maskstore4(mem, k, x);
        show(load4(mem), "store");
    }
}


void test_swap1()
{
    printf("------ test_swap1\n");
    vec4 a{ 1, 2, 3, 4};
    vec4 b{-1,-2,-3,-4};
    show(a, "a = ");
    show(b, "b = ");
    
    show(permute4(a,0x05), "permute(a,0x05)");
    show(permute4(b,0x05), "permute(b,0x05)");
    
    show(permute2f128(a,b,0x20), "permute2f128(a,b,0x20)");
    show(permute2f128(a,b,0x31), "permute2f128(a,b,0x31)");
    show(permute2f128(a,b,0x01), "permute2f128(a,b,0x01)");

    vec4 x = permute2f128(a, a, 0b00100001);
    show(shuffle4(a, x, 0b0101), "rotate >");
    show(shuffle4(x, a, 0b0101), "rotate <");
    
#if 0 //def __AVX2__
    show(rotater4(a), "rotate >");
    show(rotatel4(a), "rotate <");
#endif
    
    vec4 u = shuffle4(a, x, 0b0101);
    show(blend4(u, x, 0b1100), "rotate 3");
}


void test_swap2()
{
    printf("------ test_swap2\n");
    vec4 a{ 1,  2,  3,  4};
    vec4 b{-1, -2, -3, -4};
    show(a, "a");
    show(b, "b");
    
    show(permute2f128(a,a,0x08), "permute2f128(a,a,0x08)");
    show(permute2f128(a,a,0x21), "permute2f128(a,a,0x21)");
    show(permute2f128(a,a,0x81), "permute2f128(a,a,0x81)");
    
    show(permute2f128(a,a,0x28), "permute2f128(a,a,0x28)");
    show(permute2f128(a,a,0x81), "permute2f128(a,a,0x81)");
}


/*
 How to transform 'xyz' into 'xxx', 'yyy' and 'zzz'
 */
void test_swap4()
{
    printf("------ test_swap4\n");
    vec4 s{1, 2, 3, 4};

    {
        vec4 z = unpacklo4(s, s);
        vec4 u = unpackhi4(s, s);
        
        show(s, "src");
        show(z, "lo ");
        show(u, "hi ");
        
        show(permute2f128(z, z, 0x00), "permute2f128(z,z) 0x00");
        show(permute2f128(u, u, 0x00), "permute2f128(u,u) 0x00");
        show(permute2f128(z, z, 0x11), "permute2f128(z,z) 0x11");
        show(permute2f128(u, u, 0x11), "permute2f128(u,u) 0x11");
    }
    {
        vec4 p = permute2f128(s, s, 0x01);
        vec4 l = blend4(s, p, 0b1100);
        vec4 u = blend4(s, p, 0b0011);
        vec4 x0 = unpacklo4(l,l);
        vec4 x1 = unpackhi4(l,l);
        vec4 x2 = unpacklo4(u,u);
        vec4 x3 = unpackhi4(u,u);

        show(s, "src");
        show(x0, "xxxx");
        show(x1, "yyyy");
        show(x2, "zzzz");
        show(x3, "tttt");
    }
    {
        vec4 p = permute2f128(s, s, 0x01);
        vec4 l = blend4(s, p, 0b1100);
        vec4 u = blend4(s, p, 0b0011);
        vec4 z = unpacklo4(u,u);
        
        show(s, "xyzt");
        show(l, "xyxy");
        show(z, "zzzz");
    }
}

void test_hadd()
{
    printf("------ test_hadd\n");
    vec4 a{1, -1, 2, -2};
    vec4 b{3, -3, 4, -4};
    show(a, "a");
    show(b, "b");

    vec4 p = permute2f128(a, b, 0x20);
    vec4 q = permute2f128(a, b, 0x31);
    show(p, "p ");
    show(q, "q ");

    vec4 z = unpacklo4(p, q);
    vec4 u = unpackhi4(p, q);

    show(z, "z ");
    show(u, "u ");
}

void test_mat()
{
    printf("------ test_mat\n");
    // matrix:
    vec4 m012{0, 1, 2, -1};
    vec4 m345{3, 4, 5, -1};
    vec4 m678{6, 7, 8, -1};
    show(m012, "m012");
    show(m345, "m345");
    show(m678, "m678");
    
    // symmetrized matrix:
    vec4 z = shuffle4(m012, m345, 0b0011);
    vec4 u = permute2f128(m678, z, 0x03);
    show(z, "z");
    show(u, "u");
    
    vec4 m145 = blend4(z, m345, 0b1100);
    vec4 m258 = blend4(u, m678, 0b1100);
    show(m145, "m145");
    show(m258, "m258");
    
    // transposed matrix:
    vec4 m036 = blend4(u, shuffle4(m012, m345, 0b1000), 0b1011);
    vec4 m147 = blend4(m145, permute4(u, 0b0101), 0b0100);
    
    show(m036, "m036");
    show(m147, "m147");
}


void test_transpose2()
{
    printf("------ test_transpose2\n");
    vec4 m{0, 1, 2, 3};
    vec4 t = blend4(m, permute4(permute2f128(m,m,0x01),0b1100), 0b0110);
    show(m, "m");
    show(t, "t");
#ifdef __AVX2__
    vec4 s = permute4x64(m, 0b11011000);
    show(s, "s");
    show(permute4x64(m, 0x88), "permute4x64(m, 0x88)");
    show(permute4x64(m, 0xDD), "permute4x64(m, 0xDD)");
    show(permute4x64(m, 0x50), "permute4x64(m, 0x50)");
#endif
}


void test_transpose3()
{
    printf("------ test_transpose3\n");
    // matrix:
    vec4 m012{0, 1, 2, -6};
    vec4 m345{3, 4, 5, -1};
    vec4 m678{6, 7, 8, -1};

    show(m012, "m012");
    show(m345, "m345");
    show(m678, "m678");
    
    // symmetrized matrix:
    vec4 z = shuffle4(m012, m345, 0b0011);
    vec4 u = permute2f128(m678, z, 0x03);
    show(z, "z");
    show(u, "u");
    
    // transposed matrix:
    vec4 m036 = blend4(u, shuffle4(m012, m345, 0b1000), 0b1011);
    vec4 m258 = blend4(u, m678, 0b1100);
    vec4 m147 = blend4(z, permute4(u, 0b0101), 0b0100);

    show(m036, "m036");
    show(m147, "m147");
    show(m258, "m258");
}

void test_transpose4()
{
    printf("------ test_transpose4\n");
    // matrix:
    vec4 m0{ 0,  1,  2,  3};
    vec4 m1{ 4,  5,  6,  7};
    vec4 m2{ 8,  9, 10, 11};
    vec4 m3{12, 13, 14, 15};

    show(m0, "m0");
    show(m1, "m1");
    show(m2, "m2");
    show(m3, "m3");
    printf("\n");

    // symmetrized matrix:
    vec4 u0 = unpacklo4(m0, m1);
    vec4 u1 = unpackhi4(m1, m0);
    vec4 u2 = unpacklo4(m2, m3);
    vec4 u3 = unpackhi4(m3, m2);

    show(u0, "u0");
    show(u1, "u1");
    show(u2, "u2");
    show(u3, "u3");
    printf("\n");
    
    vec4 t0 = permute2f128(u0, u2, 0x20);
    vec4 t1 = permute2f128(u1, u3, 0x20);
    vec4 t2 = blend4(u0, permute4(u2, 0b0101), 0b0100);
    vec4 t3 = blend4(u1, permute4(u3, 0b0101), 0b0100);

    show(t0, "t0");
    show(t1, "t1");
    show(t2, "t2");
    show(t3, "t3");
}

void test_swap7()
{
    printf("------ test_swap7\n");
    vec4 s{1, 2, 3, 4};
    show(s, "source");

    show(permute4(s, 0b1010), "permute 0b1010");
    show(permute4(s, 0b0101), "permute 0b0101");

    show(shuffle4(s, s, 0b1100), "shuffle 0b1100");
    show(shuffle4(s, s, 0b0011), "shuffle 0b0011");
    
    show(permute4(s, 0b1100), "permute 0b1100");
    show(permute4(s, 0b0011), "permute 0b0011");

    show(permute4(s, 0b0000), "permute 0b0000");
    show(permute4(s, 0b1111), "permute 0b1111");

    show(unpacklo4(s, s), "unpacklo");
    show(unpackhi4(s, s), "unpackhi");

#ifdef __AVX2__
    show(permute4x64(s, 0xDD), "permute4x64 0xDD");
    show(permute4x64(s, 0x88), "permute4x64 0x88");
    show(permute4x64(s, 0xD8), "permute4x64 0xD8");
    show(permute4x64(s, 0xC9), "permute4x64 0xC9");
    show(permute4x64(s, 0xD2), "permute4x64 0xD2"); // Z X Y T
    show(permute4x64(s, 0xC9), "permute4x64 0xC9"); // Y Z X T
#endif
}

#endif

int main(int argc, char * argv[])
{
    test_swapSSE();
#ifdef __AVX__
    if ( 1 )
    {
        test_cat();
        test_deswizzle();
        test_load();
        test_broadcast();
        test_store();
    }
    if ( 1 )
    {
        //test_swap1();
        //test_swap2();
        test_swap4();
        test_hadd();
    }
    if ( 1 )
    {
        test_transpose2();
        test_transpose3();
        test_transpose4();
        test_swap7();
    }
#endif
    if ( 1 )
    {
        init();
        run(scalar,  "scalar ");
        run(vector2, "vector2");
#ifdef __AVX__
        run(vector4, "vector4");
        run(vectorU, "vectorU");
#endif
    }
    return 0;
}

