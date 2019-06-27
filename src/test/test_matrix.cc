// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/time.h>
#include <fstream>

#include "assert_macro.h"
#include "exceptions.h"
#include "random.h"
#include "tictoc.h"
#include "vecprint.h"

#include "matsparsesym.h"
#include "matsparsesym1.h"
#include "matsparsesym2.h"
#include "matsparsesymblk.h"

using namespace TicToc;

typedef MatrixSparseSymmetricBlock MatrixSparseSymmetricB;

const int N_RUN = 16;
const int N_MUL = 64;

#define PAD 4

real alpha = 2.0;
real beta = 1.0;


real checksum(int size, real const* vec, real const* ptr)
{
    real s = 0;
    for ( int kk=0; kk<size; ++kk )
        s += vec[kk] * ptr[kk];
    return s;
}


bool equal(int size, real const* a, real const* b)
{
    for ( int kk=0; kk<size; ++kk )
        if ( a != b )
            return false;
    return true;
}

// fill lower triangle of matrix
void setIndices(int fill, int*& ii, int*& jj, int mx, int bs)
{
    delete[] ii;
    delete[] jj;
    
    ii = new int[fill];
    jj = new int[fill];
    
    for ( int n = 0; n < fill; ++n )
    {
        int i = RNG.pint(mx);
        int j = RNG.pint(mx);
        ii[n] = bs * ( std::max(i,j) / bs );
        jj[n] = bs * ( std::min(i,j) / bs );
    }
}

void setVectors(size_t size, real*& x, real*& y, real*& z)
{
    free_real(x);
    free_real(y);
    free_real(z);

    x = new_real(size+PAD);
    y = new_real(size+PAD);
    z = new_real(size+PAD);

    for ( size_t n = 0; n < size; ++n )
    {
        x[n] = RNG.sreal();
        y[n] = RNG.preal();
        z[n] = 0.0;
    }
    for ( size_t n = size; n < size+PAD; ++n )
    {
        x[n] = 0.0;
        y[n] = 0.0;
        z[n] = 0.0;
    }
}

//------------------------------------------------------------------------------

template <typename MATRIXA, typename MATRIXB>
void compare(unsigned size,  MATRIXA & mat1, MATRIXB& mat2, unsigned fill)
{
    real * tmp1 = new_real(size*size);
    real * tmp2 = new_real(size*size);
    
    mat1.reset();
    mat2.reset();
    
    mat1.resize(size);
    mat2.resize(size);
    
    for ( unsigned n = 0; n < fill; ++n )
    {
        real a = 10.0 * RNG.preal();
        unsigned ii = RNG.pint(size);
        unsigned jj = RNG.pint(size);
        mat1(ii, jj) += a;
        mat2(ii, jj) += a;
    }
    
    for ( unsigned nbc = DIM; nbc < size; nbc+=DIM )
    {
        unsigned inx = DIM * ( RNG.pint(size-nbc) / DIM );
        std::clog << "Comparing matrices: size " << size << " inx " << inx << " nbc " << nbc << " ";
        
        zero_real(size*size, tmp1);
        
        mat1.addDiagonalBlock(tmp1, size, inx, nbc);
        
        //std::clog<<"mat1:\n";
        //VecPrint::print(std::clog, nbc, nbc, tmp1, size);
        
        zero_real(size*size, tmp2);
        mat2.addDiagonalBlock(tmp2, size, inx, nbc);
        
        //std::clog<<"mat2:\n";
        //VecPrint::print(std::clog, nbc, nbc, tmp2, size);
        
        real nrm = 0;
        for ( size_t n=0; n<size*size; ++n )
        {
            tmp1[n] -= tmp2[n];
            nrm += tmp1[n] * tmp1[n];
        }
        
        if ( nrm > 0 )
        {
            std::clog<<": diff:\n";
            VecPrint::print(std::clog, nbc, nbc, tmp1, size);
        }
        else
        {
            std::clog<<": identical\n";
            VecPrint::print(std::clog, nbc, nbc, tmp2, size);
        }
    }
    
    free_real(tmp1);
    free_real(tmp2);
}


template <typename MATRIX>
void fillMatrix3D(MATRIX& mat, const int ii, const int jj)
{
    mat(ii  , ii  ) += alpha;
    mat(ii+1, ii  ) -= beta;
    mat(ii+2, ii  ) -= beta;
    mat(jj  , ii  ) -= beta;
    mat(jj+1, ii  ) -= beta;
    mat(jj+2, ii  ) -= beta;
    
    mat(ii+1, ii+1) += alpha;
    mat(ii+2, ii+1) -= beta;
    mat(jj  , ii+1) -= beta;
    mat(jj+1, ii+1) -= beta;
    mat(jj+2, ii+1) -= beta;
    
    mat(ii+2, ii+2) += alpha;
    mat(jj  , ii+2) -= beta;
    mat(jj+1, ii+2) -= beta;
    mat(jj+2, ii+2) -= beta;
    
    mat(jj  , jj  ) += alpha;
    mat(jj+1, jj  ) -= beta;
    mat(jj+2, jj  ) -= beta;
    
    mat(jj+1, jj+1) += alpha;
    mat(jj+2, jj+1) -= beta;
    
    mat(jj+2, jj+2) += alpha;
}

template <typename MATRIX>
void fillMatrix2D(MATRIX& mat, const int ii, const int jj)
{
    mat(ii  , ii  ) += alpha;
    mat(ii+1, ii  ) -= beta;
    mat(jj  , ii  ) -= beta;
    mat(jj+1, ii  ) -= beta;
    
    mat(ii+1, ii+1) += alpha;
    mat(jj  , ii+1) -= beta;
    mat(jj+1, ii+1) -= beta;
    
    mat(jj  , jj  ) += alpha;
    mat(jj+1, jj  ) -= beta;
    
    mat(jj+1, jj+1) += alpha;
}

template <typename MATRIX>
void fillMatrix1D(MATRIX& mat, const int ii, const int jj)
{
    mat(ii, ii) += alpha;
    mat(jj, ii) -= beta;
    mat(jj, jj) += alpha;
}


template <typename MATRIX>
void testMatrix(MATRIX & mat,
                const int size, real const* x, real const* y, real * z,
                const int fill, int inx[], int iny[])
{
    mat.resize(size);

    tic();
    for ( int ii=0; ii<N_RUN; ++ii )
    {
        mat.reset();
        for ( int n=0; n<fill; ++n )
        {
#if ( DIM == 3 )
            fillMatrix3D(mat, iny[n], inx[n]);
#elif ( DIM == 2 )
            fillMatrix2D(mat, iny[n], inx[n]);
#else
            fillMatrix1D(mat, iny[n], inx[n]);
#endif
        }
    }
    double ts = toc();

    tic();
    for ( int ii=0; ii<N_RUN; ++ii )
    {
        mat.prepareForMultiply(1);
        for ( int n=0; n<N_MUL; ++n )
            mat.vecMulAdd(x, z);
    }
    double t1 = toc();
    
    copy_real(size, y, z);
    mat.vecMulAdd(x, z);
    real sum = checksum(size, z, x);
    
    printf("Matrix %20s : ", mat.what().c_str());
    printf("set %8.3f  mul %8.3f", ts, t1);
    printf("  checksum %+32.16f  ", sum);
    
    //---- multidimensional isotropic multiplication
    tic();
    for ( int ii=0; ii<N_RUN; ++ii )
    {
        mat.reset();
        for ( int n=0; n<fill; ++n )
            fillMatrix1D(mat, inx[n], iny[n]);
    }
    double t2 = toc();
    
    tic();
    for ( int ii=0; ii<N_RUN; ++ii )
    {
        mat.prepareForMultiply(DIM);
        for ( int n=0; n<N_MUL; ++n )
#if ( DIM >= 3 )
            mat.vecMulAddIso3D(x, z);
#else
            mat.vecMulAddIso2D(x, z);
#endif
    }
    double t3 = toc();
    
    printf("  isoset %8.3f  isomul %8.3f\n", t2, t3);
}


void testMatrices(const int size, const int fill)
{
    printf("\n **** Matrix size %i  filled %.1f %% :\n", size, fill*100.0/size/size);
    MatrixSparseSymmetric  mat0;
    MatrixSparseSymmetric1 mat1;
    MatrixSparseSymmetric2 mat2;
    MatrixSparseSymmetricB mat3;
    
    int * inx = nullptr;
    int * iny = nullptr;
    
    setIndices(fill, iny, inx, size, DIM);
    
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;

    setVectors(DIM*size, x, y, z);
    alpha = RNG.sreal();
   
    testMatrix(mat0, size, x, y, z, fill, inx, iny);
    testMatrix(mat1, size, x, y, z, fill, inx, iny);
    testMatrix(mat2, size, x, y, z, fill, inx, iny);
    testMatrix(mat3, size, x, y, z, fill, inx, iny);
    
    if ( 1 )
    {
        std::ofstream os1("mat1.txt");
        std::ofstream os3("mat3.txt");
        mat1.printSparse(os1);
        mat3.printSparse(os3);
    }

    free_real(x);
    free_real(y);
    free_real(z);
    delete[] inx;
    delete[] iny;
}


const real dir[4] = {  2, 1, -1, 3 };
const real vec[4] = { -1, 3,  1, 2 };

void fillMatrixBlock(MatrixSparseSymmetricBlock& mat, const int fill, int inx[], int iny[])
{
    SquareBlock S = SquareBlock::outerProduct(dir);
    SquareBlock U = SquareBlock::outerProduct(dir, vec);
    
    for ( int n=0; n<fill; ++n )
    {
        int ii = inx[n] - inx[n] % SquareBlock::dimension();
        int jj = iny[n] - iny[n] % SquareBlock::dimension();
        mat.diag_block(ii).sub_diag(S);
        mat.diag_block(jj).add_diag(S);
        mat.block(ii, jj).add_full(U);
    }
}

/**
This compares the Scalar and SIMD implementations of one matrix
*/
 void testMatrixBlock(MatrixSparseSymmetricBlock & mat,
                     const int size, real const* x, real const* y, real * z,
                     const int fill, int inx[], int iny[])
{
    mat.resize(size);
    
    tic();
    for ( int r = 0; r < N_RUN; ++r )
    {
        mat.reset();
        fillMatrixBlock(mat, fill, inx, iny);
    }
    double ts = toc();

    mat.prepareForMultiply(1);
    
    // calculate and compare sum for two methods:
    copy_real(size+PAD, y, z);
    mat.vecMulAdd_SCAL(x, z);
    real sum = checksum(size, z, x);
    
    copy_real(size+PAD, y, z);
    mat.vecMulAdd_SIMD(x, z);
    real res = checksum(size, z, x);

    unsigned long long time = __rdtsc();
    for ( int n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( int m=0; m<N_MUL; ++m )
            mat.vecMulAdd_SCAL(x, z);
    }
    double nop = N_MUL * N_RUN * mat.nbElements();
    double t1 = ( __rdtsc() - time ) / nop;

    time = __rdtsc();
    for ( int n = 0; n < N_RUN; ++n )
    {
        mat.prepareForMultiply(1);
        for ( int m = 0; m < N_MUL; ++m )
            mat.vecMulAdd_SIMD(x, z);
    }
    double t2 = ( __rdtsc() - time ) / nop;
    
    printf("%6i %18s ", size, mat.what().c_str());
    printf("set %8.1f scal %8.1f  simd %8.1f", ts, t1, t2);
    printf(" :  checksum  %+24.16f %+24.16f", sum, res);
    if ( sum != res )
        printf("  failed!\n");
    else
        printf("\n");
}


void testMatrixBlock(const int size, const int fill)
{
    int * inx = nullptr;
    int * iny = nullptr;
    real * x = nullptr;
    real * y = nullptr;
    real * z = nullptr;

    setIndices(fill, iny, inx, size, DIM);
    setVectors(size, x, y, z);
    alpha = RNG.sreal();
    
    MatrixSparseSymmetricBlock mat;
    testMatrixBlock(mat, size, x, y, z, fill, inx, iny);
    
    free_real(x);
    free_real(y);
    free_real(z);
    delete[] inx;
    delete[] iny;
}


 int main( int argc, char* argv[] )
{
    RNG.seed();
    if ( 0 )
    {
        // small tests to check correctness:
        MatrixSparseSymmetric1 mat1;
        MatrixSparseSymmetricB mat3;
        compare(4*3, mat1, mat3, 1<<4);
        compare(4*7, mat1, mat3, 1<<5);
        compare(4*11, mat1, mat3, 1<<6);
        compare(4*33, mat1, mat3, 1<<16);
    }
    if ( 0 )
    {
        testMatrices(6, 1);
        testMatrices(12, 1);
        testMatrices(DIM*7, 1);
        testMatrices(DIM*17, 2);
        testMatrices(DIM*33, 1111);
    }
    if ( 0 )
    {
        printf("\ntest_matrix BLOCK_SIZE %i  (BLD %i)\n", SquareBlock::dimension(), SquareBlock::stride());
        size_t siz = DIM;
        for ( int i = 0; i < 14; ++i )
        {
            siz = 1 + ( 2 << i ) * RNG.preal();
            size_t fill = 1 + siz * siz * 0.01 * RNG.preal();
            testMatrixBlock(siz, fill);
        }
    }
    if ( 0 )
    {
        testMatrixBlock(DIM*17, 2);
        testMatrixBlock(DIM*347, 1019);
        testMatrixBlock(DIM*753, 43039);
    }
    if ( 0 )
    {
        testMatrixBlock(2253, 1<<14);
        testMatrixBlock(2253, 1<<14);
        testMatrixBlock(2253, 1<<14);
        testMatrixBlock(2253, 1<<14);
    }
    if ( 0 )
    {
        testMatrixBlock(DIM*1251, 25821);
        testMatrixBlock(DIM*1785, 153034);
        testMatrixBlock(DIM*2311, 231111);
        //testMatrixBlock(DIM*3217, 671234);
    }
    if ( 1 )
    {
        //testMatrices(DIM*17, 23);
        testMatrices(DIM*91, 1<<12);
        testMatrices(DIM*197, 1<<14);
        testMatrices(DIM*437, 1<<17);
        testMatrices(DIM*713, 1<<18);
    }
    return EXIT_SUCCESS;
}


