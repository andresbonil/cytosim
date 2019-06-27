// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

/**
 * ------------------------------------------------------------------------------
 *                   -- Meca is the heart of Cytosim --
 * ------------------------------------------------------------------------------
 *             It solves the equations of motion for the Mecables,
 *      using implicit integration and iterative methods with sparse matrix
 * ------------------------------------------------------------------------------
 * @todo See if Lagrangian dynamics could work better than constrainted dynamics
 * @todo Implement the PARDISO sparse matrix format
 * @todo Check if IDR(s) would perform better than BCGS or GMRES
 * ------------------------------------------------------------------------------
 */

#include <fstream>

#include "meca.h"
#include "mecable.h"
#include "messages.h"
#include "simul_prop.h"
#include "cblas.h"
#include "clapack.h"
#include "exceptions.h"
#include "vecprint.h"
#include "filepath.h"
#include "tictoc.h"
#include "bicgstab.h"
#include "gmres.h"

#include "meca_inter.cc"

/**
 Add correction term to the constrainted dynamics
 The effect is to stabilize fibers under traction, at some modest CPU cost.
*/
#define ADD_PROJECTION_DIFF 1


/**
 The forces are usually:
 
     force = vBAS + ( mB + mC + mR + mdiffP ) * vPTS
 
 where mR represent the bending elasticity of fibers.
*/

/// this define will enable explicit integration (should be off)
#define EXPLICIT_INTEGRATION 0


/// this is used to check the validity of the results
#define not_a_number(x) ((x) != (x))


/// use this to generate code for validation
#define DEBUG_MECA 0


/// number of threads running in parallel
#define NUM_THREADS 1


#if NUM_THREADS > 1
// Parallelization uses Intel's OpenMP
#include <omp.h>
#endif


//------------------------------------------------------------------------------
#pragma mark - Allocate

Meca::Meca()
: objs(32)
{
    ready_ = 0;
    nbPts = 0;
    allocated_ = 0;
    vPTS = nullptr;
    vSOL = nullptr;
    vBAS = nullptr;
    vRND = nullptr;
    vRHS = nullptr;
    vFOR = nullptr;
    vTMP = nullptr;
    vMEM = nullptr;
    useMatrixC = false;
    drawLinks = false;
}


void allocate_vector(size_t s, real *& ptr, bool reset)
{
    free_real(ptr);
    ptr = new_real(s);
    if ( reset )
        zero_real(s, ptr);
}


void Meca::allocate(size_t alc)
{
    //allocate the vectors
    if ( alc > allocated_ )
    {
        // make a multiple of chunk to align pointers:
        allocated_ = chunk_real(alc);
        
        // pad with 4 doubles to allow SIMD instruction burr
        alc = DIM * allocated_ + 4;
        allocate_vector(alc, vPTS, 1);
        allocate_vector(alc, vSOL, 1);
        allocate_vector(alc, vBAS, 0);
        allocate_vector(alc, vRND, 1);
        allocate_vector(alc, vRHS, 1);
        allocate_vector(alc, vFOR, 1);
        allocate_vector(alc, vTMP, 0);
#if NUM_THREADS > 1
        allocate_vector(alc, vMEM, 0);
#endif
    }
}


void Meca::release()
{
    //std::clog << "Meca::release()\n";
    free_real(vPTS);
    free_real(vSOL);
    free_real(vBAS);
    free_real(vRND);
    free_real(vRHS);
    free_real(vFOR);
    free_real(vTMP);
    free_real(vMEM);
    vPTS = nullptr;
    vSOL = nullptr;
    vBAS = nullptr;
    vRND = nullptr;
    vRHS = nullptr;
    vFOR = nullptr;
    vTMP = nullptr;
    vMEM = nullptr;
}


unsigned Meca::largestMecable() const
{
    unsigned res = 0;
    for ( Mecable * mec : objs )
        res = std::max(res, mec->nbPoints());
    return res;
}

//------------------------------------------------------------------------------
#pragma mark - Multiply


// shortcut

#if ( DIM == 1 )
#   define VECMULADDISO  vecMulAdd
#elif ( DIM == 2 )
#   define VECMULADDISO  vecMulAddIso2D
#elif ( DIM == 3 )
#   define VECMULADDISO  vecMulAddIso3D
#endif


/**
 calculate the forces into `F`, given the Mecable coordinates `X`:
 
     F <- B + mB * X + mC * X

 If B == 0, this term is ommited. With B = vBAS and X = vPTS, the procedure
 calculates the forces in the system in `F`:
 
     F <- vBAS + mB * vPTS + mC * vPTS

 */
#if NUM_THREADS > 1
void Meca::calculateForces(const real* X, real const* B, real* F) const
{
    assert_true( empty() || ( X != F && X != B && F != B ));

    /*
    We partition the matrix into equal parts, assuming it is of constant density,
    symmetric and stored in its lower triangle, such that columns of low index
    have more elements: n_elements(inx) ~ ( size - inx )
    Hence to split in two we use these ratios:
        s = 1.0 - sqrt(1/2);
    and to split in 3 equal parts:
        s1 = 1.0 - sqrt(2/3);
        s2 = 1.0 - sqrt(1/3);
    or to split in 4 equal parts:
        s1 = 1.0 - sqrt(3/4);
        s2 = 1.0 - sqrt(2/4);
        s3 = 1.0 - sqrt(1/4);
     */
    index_t dim = dimension();

    if ( useMatrixC )
    {
        index_t spl1 = index_t( dim * 0.184 );
        index_t spl2 = index_t( dim * 0.423 );
#if ( 0 )
        std::clog << " matrixC split [ 0 " << spl1 << " " << spl2 << " " << dim <<" ] :  ";
        std::clog << mC.nbElements(0, spl1) << " " << mC.nbElements(spl1, spl2) << " " << mC.nbElements(spl2, siz) << "\n";
#endif
        #pragma omp parallel sections num_threads(4)
        {
            #pragma omp section
            {
                // F <- B + mB * X
                if ( B )
                    copy_real(dim, B, F);      //blas::xcopy(dimension(), B, 1, F, 1);
                else
                    zero_real(dim, F);
                
                mB.VECMULADDISO(X, F);
            }
            #pragma omp section
            {
                zero_real(dim, vRND);
                mC.vecMulAdd(X, vRND, 0, spl1);
            }
            #pragma omp section
            {
                zero_real(dim, vTMP);
                mC.vecMulAdd(X, vTMP, spl1, spl2);
            }
            #pragma omp section
            {
                zero_real(dim, vMEM);
                mC.vecMulAdd(X, vMEM, spl2, dim);
            }
        }
    }
    else
    {
        index_t nbp = nbPts;
        index_t spl1 = index_t( nbp * 0.134 );
        index_t spl2 = index_t( nbp * 0.293 );
        index_t spl3 = index_t( nbp * 0.5 );
#if ( 0 )
        std::clog << " matrixB split [ 0 " << spl1 << " " << spl2 << " " << spl3 << " " << nbp <<" ] :  ";
        std::clog << mB.nbElements(0, spl1) << " " << mB.nbElements(spl1, spl2) << " ";
        std::clog << mB.nbElements(spl2, spl3) << " " << mB.nbElements(spl3, siz) << "\n";
#endif
        #pragma omp parallel sections num_threads(4)
        {
            #pragma omp section
            {
                // F <- B + mB * X
                if ( B )
                    copy_real(dim, B, F);
                else
                    zero_real(dim, F);
                mB.VECMULADDISO(X, F, 0, spl1);
            }
            #pragma omp section
            {
                zero_real(dim, vRND);
                mB.VECMULADDISO(X, vRND, spl1, spl2);
            }
            #pragma omp section
            {
                zero_real(dim, vTMP);
                mB.VECMULADDISO(X, vTMP, spl2, spl3);
            }
            #pragma omp section
            {
                zero_real(dim, vMEM);
                mB.VECMULADDISO(X, vMEM, spl3, nbp);
            }
        }
    }

    // finally sum up all contributions into F
    #pragma vector aligned
    for ( index_t i = 0; i < dim; ++i )
        F[i] = ( F[i] + vRND[i] ) + ( vTMP[i] + vMEM[i] );
}

#else

/**
 calculate the forces into `F`, given the Mecable coordinates `X`:
 
     F <- B + mB * X + mC * X

 If `B == 0`, this term is ommited. With `B = vBAS` and `X = vPTS`, this
 effectively calculates the forces in the system in `F`:
 
     F <- vBAS + mB * vPTS + mC * vPTS

 */
void Meca::calculateForces(const real* X, real const* B, real* F) const
{
    assert_true( empty() || ( X != F && X != B && F != B ));

    // F <- B  of F <- 0
    if ( B )
        copy_real(dimension(), B, F);      //blas::xcopy(dimension(), B, 1, F, 1);
    else
        zero_real(dimension(), F);
    
    // F <- F + mB * X
#if ( DIM == 1 )
    mB.vecMulAdd(X, F);
#elif ( DIM == 2 )
    mB.vecMulAddIso2D(X, F);
#elif ( DIM == 3 )
    mB.vecMulAddIso3D(X, F);
#endif

    if ( useMatrixC )
    {
        // F <- F + mC * X
        mC.vecMulAdd(X, F);
    }
}
#endif


void Meca::addAllRigidity(const real* X, real* Y) const
{
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            const index_t inx = DIM * (*mci)->matIndex();
            (*mci)->addRigidity(X+inx, Y+inx);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : objs )
    {
        const index_t inx = DIM * mec->matIndex();
        mec->addRigidity(X+inx, Y+inx);
    }
#endif
}


/**
 calculate the matrix vector product corresponding to 'mec'
 
     Y <- X + alpha * speed( Y + P' * X );
 
 */
inline void multiply1(Mecable const* mec, real alpha, const real* xxx, real* yyy)
{
#if ( DIM > 1 )
    mec->addRigidity(xxx, yyy);
#endif

#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(xxx, yyy);
#endif

    mec->setSpeedsFromForces(yyy, alpha, yyy);
    
    /*
     Bypass the Projection (NEVER ENABLE THIS!):
     real alpha = mec->nbPoints() * time_step / mec->dragCoefficient();
     blas::xaxpy(bs, -alpha, tmp, 1, yyy, 1);
     PRINT_ONCE("CRAZY: Projection is bypassed!!!");
     */

    // Y <- Y + X
    //blas::xaxpy(DIM * mec->nbPoints(), 1.0, xxx, 1, yyy, 1);
    blas::add(DIM*mec->nbPoints(), xxx, yyy);
}


/**
 calculate the matrix product needed for the conjugate gradient algorithm
 
     Y <- X - time_step * speed( mB + mC + P' ) * X;
 
 */
void Meca::multiply( const real* X, real* Y ) const
{
    // Y <- ( mB + mC ) * X
    calculateForces(X, nullptr, Y);
    
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            const index_t inx = DIM * (*mci)->matIndex();
            multiply1(*mci, -time_step, X+inx, Y+inx);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : objs )
    {
        const index_t inx = DIM * mec->matIndex();
        multiply1(mec, -time_step, X+inx, Y+inx);
    }
#endif
}


/**
 calculate the matrix vector product corresponding to 'mec'
 
     X <- PRECOND * [ X + alpha * speed( Y + meca_forces ) ]
 
 */
inline void multiply1P(Mecable const* mec, real alpha, real* xxx, real* yyy)
{
    const int bks = DIM * mec->nbPoints();

#if ( DIM > 1 )
    mec->addRigidity(xxx, yyy);
#endif

#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
        mec->addProjectionDiff(xxx, yyy);
#endif

    mec->setSpeedsFromForces(yyy, alpha, yyy);
    
    // X <- Y + X
    //blas::xaxpy(bks, 1.0, yyy, 1, xxx, 1);
    blas::add(bks, yyy, xxx);
    
    if ( mec->useBlock() )
    {
        int info = 0;
        lapack::xgetrs('N', bks, 1, mec->block(), bks, mec->pivot(), xxx, bks, &info);
    }
}

/**
 This is similar to
 
     multiply(X, T);       // T <- M*X
     precondition(T, X);   // X <- P*T
 
 if the value of 'T' is not used
 */
void Meca::multiplyP(real* X, real* T) const
{
    // T <- ( mB + mC ) * X
    calculateForces(X, nullptr, T);
    
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            const index_t inx = DIM * (*mci)->matIndex();
            multiply1P(*mci, -time_step, X+inx, T+inx);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : objs )
    {
        const index_t inx = DIM * mec->matIndex();
        multiply1P(mec, -time_step, X+inx, T+inx);
    }
#endif
}

//------------------------------------------------------------------------------
#pragma mark - Helper functions


/**
 Fill-in matrix 'dst' as the duplicate of 'src', for each 'DIM' dimension.
 For 'DIM==1', this makes a simple copy
 For 'DIM==2', 'src' is copied twice, into odd indices, and into even indices.
 For 'DIM==3', three copies of 'src' are made into 'dst'.
 
 Both 'src' and 'dst' must be symmetrix square matrices. 
 The size of 'dst' is DIM times the size of 'src'.
 Only the upper diagonal of 'src' is specified.
 Matrix Y is specified in full.
 */

void duplicate_matrix(unsigned siz, real const* src, real * dst)
{
    unsigned ddd = DIM * siz;
    
    zero_real(ddd*ddd, dst);
    
    for ( unsigned ii = 0; ii < siz; ++ii )
    {
        real xx = src[ii + siz * ii];
        
        unsigned kk = ( ddd+1 ) * DIM * ii;
        for ( unsigned d = 0; d < DIM; ++d, kk += ddd+1 )
            dst[kk] = xx;
        
        for ( unsigned jj = ii+1; jj < siz; ++jj )
        {
            xx = src[ii + siz * jj];
            kk = DIM * ( ii + ddd * jj );
            unsigned ll = DIM * ( jj + ddd * ii );
            for ( unsigned d = 0; d < DIM; ++d, kk += ddd+1, ll += ddd+1 )
            {
                dst[kk] = xx;
                dst[ll] = xx;
            }
        }
    }
    
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, src);
    std::clog << "Duplicated:\n";
    VecPrint::print(std::clog, ddd, ddd, dst);
#endif
}


/**
 This should symmetrize a matrix, and also copy the terms
 that are within the first subspace `X` into the other dimensions
 input: upper triangular matrix
 ouput: full symmetric matrix
 */
void expand_matrix(unsigned siz, real * mat)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
    
    for ( unsigned jj = 0; jj < siz; jj += DIM  )
    {
        for ( unsigned ii = 0; ii < jj; ii += DIM  )
        {
            real val = mat[ii+siz*jj];
            // expand term in other dimensions:
            for ( int d = 1; d < DIM; ++d )
                mat[ii+d+siz*(jj+d)] = val;
            
            // symmetrize matrix:
            for ( int d = 0; d < DIM; ++d )
                mat[jj+d+siz*(ii+d)] = val;
        }
        // expand diagonal term in other dimensions:
        real val = mat[jj+siz*jj];
        for ( int d = 1; d < DIM; ++d )
            mat[jj+d+siz*(jj+d)] = val;
    }

#if ( 0 )
    std::clog << "Expanded:\n";
    VecPrint::print(std::clog, siz, siz, mat, siz);
#endif
}


/**
 Set to zero all the terms that are not within 'diag' from the diagonal.
 With 'diag==0' the entire matrix is set to zero.
 With 'diag==1', only the diagonal is kept.
 With 'diag==2', the matrix is made tri-diagonal
 etc.
 */
void truncate_matrix(unsigned siz, real* mat, unsigned diag)
{
#if ( 0 )
    std::clog << "\nOriginal:\n";
    VecPrint::print(std::clog, siz, siz, mat);
#endif
    
    for ( unsigned ii = 0; ii < siz; ++ii )
    {
        real * col = mat + siz * ii;
        for ( unsigned jj = 0; jj+diag < ii; ++jj )
            col[jj] = 0.0;
        
        for ( unsigned kk = ii+diag+1; kk < siz; ++kk )
            col[kk] = 0.0;
    }
    
#if ( 0 )
    std::clog << "Truncated:\n";
    VecPrint::print(std::clog, siz, siz, mat);
#endif
}


/// sim(element^2) / sum(diagonal^2)
real off_diagonal_norm(int siz, real * mat)
{
    real all = 0;
    for ( int k = 0; k < siz*siz; ++k )
        all += mat[k] * mat[k];

    real dia = 0;
    for ( int k = 0; k < siz*siz; k+=siz+1 )
        dia += mat[k] * mat[k];

    return sqrt( ( all - dia ) / dia );
}


/// set all values between '-val' and 'val' to zero
void threshold_matrix(int siz, real * mat, real val)
{
    for ( int k = 0; k < siz*siz; ++k )
    {
        if ( fabs(mat[k]) < val )
            mat[k] = 0.0;
    }
}


/// set 'mat' of order `siz` to `diag * I`
void diagonal_matrix(int siz, real * mat, real val)
{
    for ( int k = 0; k < siz*siz; ++k )
        mat[k] = 0.0;
    for ( int k = 0; k < siz*siz; k+=siz+1 )
        mat[k] = val;
}


/// erase all off-diagonal terms in `mat` of order `siz`
void make_diagonal(int siz, real * mat)
{
    for ( int j = 0; j < siz; ++j )
    {
        real * col = mat + j * siz;
        int i;
        for ( i = 0; i < j; ++i )
            col[i] = 0.0;
        for ( i = j+1; i < siz; ++i )
            col[i] = 0.0;
    }
}


/// a test matrix with integer components
void test_matrix(int siz, real * mat)
{
    for ( int i = 0; i < siz; ++i )
    for ( int j = 0; j < siz; ++j )
        mat[i+siz*j] = j - i;
}


/**
 Convert a full matrix into a LAPACK banded matrix data suitable for LU factorization.
 `src` is a square matrix of side 'siz'
 `dst` is of size ldd * siz, with ldd > ku+2*kl+1
 
 with indices starting at 1 (LAPACK documentation):
 src(i,j) is stored in dst(ku+1+i-j, j) for max(1, j-ku) <= i <= min(siz, j+kl)

 considering that ku <- ku + kl:
 src(i,j) is stored in dst(ku+kl+1+i-j, j) for max(1, j-ku) <= i <= min(siz, j+kl)

 with indices starting at 0:
 src(i,j) is stored in dst(ku+kl+i-j, j) for max(0, j-ku) <= i <= min(siz-1, j+kl)

 */
void banded_matrix(int siz, real const* src, int kl, int ku, real * dst, int ldd)
{
    assert_true( ldd == kl+kl+ku+1 );
#if ( 0 )
    if ( siz < 64 )
    {
        std::clog << "\noriginal:\n";
        VecPrint::print(std::clog, siz, siz, src, siz);
    }
#endif
    for ( int jj = 0; jj < siz; ++jj )
    {
        // dst[ii+ldd-kl-1-jj+ldd*jj] = src[ii+siz*jj];
        int add = ldd - kl - 1 - jj + ldd * jj;
        int inf = add + std::max(0, jj-ku);
        int sup = add + std::min(siz-1, jj+kl);
        int off = siz * jj - add;

        //std::clog << " inf " << inf << " sup " << sup << " off " << off << " jj " << jj << "  shift " << shift << "\n";
        //std::clog << " dst[ " << inf << " " << sup << " ] <- src[ " << inf+off << " " << sup+off << " ] \n";
        
        if ( off > 0 )
        {
            for (int ii = inf; ii <= sup; ++ii )
                dst[ii] = src[ii+off];
        }
        else
        {
            for (int ii = sup; ii >= inf; --ii )
                dst[ii] = src[ii+off];
        }
    }
    // zero out values:
    for ( int jj = 0; jj < siz; ++jj )
    {
        int add = ldd - kl - 1 - jj;
        int inf = add + std::max(0, jj-ku);
        int sup = add + std::min(siz-1, jj+kl);

        real * col = dst + ldd * jj;
        for (int ii = 0; ii < inf; ++ii )
            col[ii] = 0.0;
        for (int ii = sup+1; ii < ldd; ++ii )
            col[ii] = 0.0;
    }
#if ( 0 )
    if ( siz < 64 )
    {
        std::clog << " banded_matrix (size " << siz << ") :\n";
        VecPrint::print(std::clog, ldd, siz, dst, ldd);
    }
#endif
}


/*
 uses power iterations to estimate the largest eigenvalue of `mat * tam + alpha * I`
 Vector `vec` is used to initialize the algorithm
 @returns an estimate of the largest eigenvalue
 The precision of the estimate is low: 10%
 */
real largest_eigenvalue(int siz, real const* blk, int const* piv, real const* mat, real alpha, real * vec, real * tmp)
{
    assert_true(siz > 0);
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    //fprintf(stderr, "      power size %i eig %10.6f\n", siz, eig);

    int n, info = 0;
    for ( n = 0; n < siz; n += 2 )
    {
        blas::xcopy(siz, vec, 1, tmp, 1);
        lapack::xgetrs('N', siz, 1, blk, siz, piv, tmp, siz, &info);
        assert_true(info==0);
        blas::xgemv('N', siz, siz, 1.0/eig, mat, siz, tmp, 1, alpha/eig, vec, 1);
        oge = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16, siz), vec, 3);
        
        blas::xcopy(siz, vec, 1, tmp, 1);
        lapack::xgetrs('N', siz, 1, blk, siz, piv, tmp, siz, &info);
        assert_true(info==0);
        blas::xgemv('N', siz, siz, 1.0/oge, mat, siz, tmp, 1, alpha/oge, vec, 1);
        eig = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16, siz), vec, 3);
        
        //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
        
        if ( fabs(oge-eig) < TOLERANCE * ( fabs(eig) + fabs(oge) ) )
            break;
    }
    //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
    
    return std::max(eig, oge);
}


/*
 uses power iterations to estimate the largest eigenvalue of `mat * tam + alpha * I`
 Vector `vec` is used to initialize the algorithm
 @returns an estimate of the largest eigenvalue
 The precision of the estimate is low: 10%
 */
real largest_eigenvalue(int siz, real const* mat, real const* tam, real alpha, real * vec, real * tmp)
{
    const real TOLERANCE = 0.05;
    real oge, eig = blas::nrm2(siz, vec);
    
    int n;
    for ( n = 0; n < siz; n += 2 )
    {
        blas::xgemv('N', siz, siz, 1.0/eig, mat, siz, vec, 1,       0.0, tmp, 1);
        blas::xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, alpha/eig, vec, 1);
        oge = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16, siz), vec, 3);
        
        blas::xgemv('N', siz, siz, 1.0/oge, mat, siz, vec, 1,       0.0, tmp, 1);
        blas::xgemv('N', siz, siz, 1.0,     tam, siz, tmp, 1, alpha/oge, vec, 1);
        eig = blas::nrm2(siz, vec);
        //VecPrint::print(std::clog, std::min(16, siz), vec, 3);
        //fprintf(stderr, "      power iter %3i: eigen %10.6f %10.6f\n", n, eig, oge);
        if ( fabs(oge-eig) < TOLERANCE * ( fabs(eig) + fabs(oge) ) )
            break;
    }
    //fprintf(stderr, "      power size %4i iter %3i: eigen %10.6f %10.6f\n", siz, n, eig, oge);
    
    return std::max(eig, oge);
}

//------------------------------------------------------------------------------
#pragma mark - Precondition


/**
 Return the total diagonal block corresponding to an Object, which is:
 
     I - time_step * P ( mB + mC + P' )
 
 This block is square but not symmetric!
 */
void Meca::extractBlock(real* res, const Mecable * mec) const
{
    const unsigned ps = mec->nbPoints();
    const unsigned bs = DIM * ps;
    
    zero_real(bs*bs, res);
    
#if ( DIM > 1 )
    // set the Rigidity terms:
    mec->addRigidityUpper(res, bs);
    //std::clog<<"Rigidity block " << mec->reference() << "\n";
    //VecPrint::print(std::clog, bs, bs, tmp, bs);
#endif
    
    mB.addTriangularBlock(res, bs, mec->matIndex(), ps, DIM);
    
    expand_matrix(bs, res);
    
    if ( useMatrixC )
        mC.addDiagonalBlock(res, bs, DIM*mec->matIndex(), bs);
    
#if ( 0 )
    std::clog<<"mB+mC block:\n";
    VecPrint::print(std::clog, bs, bs, res, bs);
#endif
    
#if ADD_PROJECTION_DIFF
    if ( mec->hasProjectionDiff() )
    {
        // Include the corrections P' in preconditioner, vector by vector.
        real* tmp = vTMP + DIM * mec->matIndex();
        zero_real(bs, tmp);
        for ( unsigned ii = 0; ii < bs; ++ii )
        {
            tmp[ii] = 1.0;
            mec->addProjectionDiff(tmp, res+bs*ii);
            tmp[ii] = 0.0;
        }
        //std::clog<<"dynamic with P'\n";
        //VecPrint::print(std::clog, bs, bs, res, bs);
    }
#endif
    
    //compute the projection, by applying it to each column vector:
    /*
     This could be vectorized by having setSpeedsFromForces()
     accept multiple vectors as arguments, using SIMD instructions
     */
    for ( unsigned ii = 0; ii < bs; ++ii )
        mec->setSpeedsFromForces(res+bs*ii, -time_step, res+bs*ii);
    
    // add Identity matrix:
    for ( unsigned ii = 0; ii < bs; ++ii )
        res[bs*ii+ii] += 1.0;
}


/**
 This version builds the diagonal block directly from Meca:multiply().
 This is a very slow method that calls 'multiply()' n-times, where
 'n' is the size of the block. It does not use any of the Meca vectors.
 
 It should be used for validation only.
*/
void Meca::extractBlockSlow(real* res, const Mecable* mec) const
{
    const unsigned dim = dimension();
    const unsigned bks = DIM * mec->nbPoints();
    const unsigned off = DIM * mec->matIndex();
    
    assert_true( off+bks <= dim );
    real * vec = new_real(dim);
    real * tmp = new_real(dim);
    
    zero_real(dim, vec);
    //zero_real(bs*bs, res);
    
    // proceed column by column:
    for ( unsigned jj = 0; jj < bks; ++jj )
    {
        vec[jj+off] = 1.0;
        multiply(vec, tmp);
        vec[jj+off] = 0.0;
        blas::xcopy(bks, tmp+off, 1, res+jj*bks, 1);
    }
    
    free_real(vec);
    free_real(tmp);
}


// DEBUG: compare `blk` with block extracted using extractBlockSlow()
void Meca::compareBlocks(const Mecable * mec, const real* blk)
{
    const unsigned bs = DIM * mec->nbPoints();
    real * wrk = new_real(bs*bs);
    
    extractBlockSlow(wrk, mec);
    
    blas::xaxpy(bs*bs, -1.0, blk, 1, wrk, 1);
    real err = blas::nrm2(bs*bs, wrk);
 
    std::clog << "compareBlocks ";
    std::clog << std::setw(10) << mec->reference() << " " << std::setw(6) << bs;
    std::clog << "  | B - K | = " << std::setprecision(3) << err << std::endl;
    
    if ( err > 100 * REAL_EPSILON )
    {
        VecPrint::sparse(std::clog, bs, bs, wrk, bs, 3, (real)0.1);
        
        int s = std::min(16U, bs);
        extractBlockSlow(wrk, mec);
        std::clog << " blockS\n";
        VecPrint::print(std::clog, s, s, wrk, bs, 3);
        
        std::clog << " block \n";
        VecPrint::print(std::clog, s, s, blk, bs, 3);
    }
    
    free_real(wrk);
}


/**
 Multiply here `blk` with the dynamic block extracted by extractBlockSlow()
 and check that we recover the identity matrix
 */
void Meca::testBlock(const Mecable * mec, const real* blk)
{
    const unsigned bs = DIM * mec->nbPoints();
    
    std::clog << "  testBlock " << mec->useBlock() << " ";
    std::clog << std::setw(10) << mec->reference() << " " << std::setw(6) << bs;
 
    if ( !mec->useBlock() )
    {
        std::clog << std::endl;
        return;
    }
    
    real * wrk = new_real(bs*bs);
    real * tmp = new_real(bs*bs);
    real * vec = new_real(bs);
    
    extractBlockSlow(wrk, mec);
   
    int info = 0;
    blas::xcopy(bs*bs, wrk, 1, tmp, 1);
    if ( mec->useBlock() )
        lapack::xgetrs('N', bs, bs, blk, bs, mec->pivot(), tmp, bs, &info);
    
    for ( unsigned k=0; k < bs*bs; k += 1+bs )
        tmp[k] -= 1.0;
    
    real err = blas::nrm2(bs*bs,tmp) / bs;
    std::clog << " | 1 - PM | = " << std::setprecision(3) << err;
    
    if ( 1 )
    {
        // chose initial vector for power iteration
        blas::xcopy(bs, vRHS+DIM*mec->matIndex(), 1, vec, 1);
        //this is not valid for triangular preconditionner
        real eig = -1;
        if ( mec->useBlock() )
            eig = largest_eigenvalue(bs, blk, mec->pivot(), wrk, -1.0, vec, tmp);
        std::clog << "  eigen(PM) = " << eig;
    }
    
    std::clog << std::endl;

    if ( err > 1 )
    {
        // print preconditionner block for visual inspection:
        unsigned s = std::min(16U, bs);
        std::clog << " matrix: \n";
        VecPrint::print(std::clog, s, s, wrk, bs);
        std::clog << "\nprecond: \n";
        VecPrint::print(std::clog, s, s, blk, bs);
        std::clog << "\nprecond * matrix:\n";
        VecPrint::print(std::clog, s, s, tmp, bs);
        std::clog << "\n";
    }
    free_real(vec);
    free_real(tmp);
    free_real(wrk);
}


/**
Compute block of the preconditionner corresponding to 'mec'
 */
void Meca::computePreconditionner(Mecable* mec)
{
    unsigned bs = DIM * mec->nbPoints();

    mec->allocateBlock();
    real* blk = mec->block();
 
    // extract diagonal matrix block corresponding to this Mecable:
    extractBlock(blk, mec);

    //compareBlocks(mec, blk);

    // calculate LU factorization:
    int info = 0;
    lapack::xgetf2(bs, bs, blk, bs, mec->pivot(), &info);
    
    if ( info == 0 )
    {
        mec->useBlock(1);
        //testBlock(mec, blk);
        //std::clog << "Meca::computePreconditionner(" << mec->reference() << ")\n";
    }
    else
    {
        assert_true(mec->useBlock() == 0);
        std::clog << "Meca::computePreconditionner failed (lapack::xgetf2, info " << info << ")\n";
    }
}


/// Compute all the blocks of the preconditionner
/**
 This is method 1, that should perform well and can be multithreaded
 */
void Meca::computePreconditionner()
{
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            computePreconditionner(*mci);
            mci += NUM_THREADS;
        }
        //printf("thread %i complete %i\n", omp_get_thread_num(), TicToc::microseconds());
    }
#else
    for ( Mecable * mec : objs )
        computePreconditionner(mec);
#endif
}


void Meca::precondition(const real* X, real* Y) const
{
    assert_true( X != Y );
    
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        int info;
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            Mecable const* mec = *mci;
            const int bs = DIM * mec->nbPoints();
            real const* xxx = X + DIM * mec->matIndex();
            real * yyy = Y + DIM * mec->matIndex();
            blas::xcopy(bs, xxx, 1, yyy, 1);
            if ( mec->useBlock() )
                lapack::xgetrs('N', bs, 1, mec->block(), bs, mec->pivot(), yyy, bs, &info);
            mci += NUM_THREADS;
        }
    }
#else
    blas::xcopy(dimension(), X, 1, Y, 1);
    for ( Mecable const* mec : objs )
    {
        const int bs  = DIM * mec->nbPoints();
        const int inx = DIM * mec->matIndex();
        int info = 0;
        if ( mec->useBlock() )
            lapack::xgetrs('N', bs, 1, mec->block(), bs, mec->pivot(), Y+inx, bs, &info);
    }
#endif
}


//------------------------------------------------------------------------------
#pragma mark - Solve


/// function to sort Mecables
int smaller_mecable(const void * ap, const void * bp)
{
    Mecable const** a = (Mecable const**)(ap);
    Mecable const** b = (Mecable const**)(bp);

    if ( (*a)->nbPoints() > (*b)->nbPoints() ) return -1;
    if ( (*a)->nbPoints() < (*b)->nbPoints() ) return  1;
    return 0;
}


/**
 Allocate and reset matrices and vectors necessary for Meca::solve(),
 copy coordinates of Mecables into vPTS[]
 */
void Meca::prepare(SimulProp const* prop)
{
#if NUM_THREADS > 1
    /*
     Sorting Mecables can improve multithreaded performance by distributing
     the work more equally between threads. Note that his operation is not free
     and for large systems random partitionning may not be so bad. Moreover for
     homogeneous systems (if all filaments have the same length) this is useless.
    */
    objs.sort(smaller_mecable);
    
    /*
    for ( Mecable const* mec : objs )
        std::clog << mec->reference() << " " << mec->nbPoints() << "\n";
     */
#endif
    
    /*
     Attributes the position in the vector/matrix to each Mecable
     */
    index_t cnt = 0;
    for ( Mecable * mec : objs )
    {
        mec->matIndex(cnt);
        cnt += mec->nbPoints();
    }
    nbPts = cnt;
    allocate(cnt);
    
    //allocate the sparse matrices:
    mB.resize(cnt);
    mB.reset();

    mC.resize(DIM*cnt);
    mC.reset();
    
    // reset base:
    zero_real(DIM*cnt, vBAS);
    
    // get global time step
    time_step = prop->time_step;
    
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            Mecable * mec = *mci;
            mec->putPoints(vPTS+DIM*mec->matIndex());
            mec->prepareMecable();
            mec->useBlock(0);
            mci += NUM_THREADS;
        }
    }
#else
    for ( Mecable * mec : objs )
    {
        mec->putPoints(vPTS+DIM*mec->matIndex());
        mec->prepareMecable();
        mec->useBlock(0);
    }
#endif
}


/**
 Prepare matrices mB and mC for multiplication
 This should be called after setInteractions()
 */
void Meca::prepareMatrices()
{
    mB.prepareForMultiply(DIM);
    
    if ( mC.nonZero() )
    {
        useMatrixC = true;
        mC.prepareForMultiply(1);
    }
    else
        useMatrixC = false;
}


/**
 Calculates the force in the objects, that can be accessed by Mecable::netForce()
 and calculate the speed of the objects in vRHS, in the abscence of Thermal motion,
 ie. the motion is purely due to external forces.

 This also sets the Lagrange multipliers for the Fiber.
 
 The function will not change the position of the Mecables.
 */
void Meca::computeForces()
{
    prepareMatrices();
    
    // calculate forces in vFOR, but without adding Rigidity and Brownian noise:
    calculateForces(vPTS, vBAS, vFOR);
    
    // return forces to Mecable, and compute Lagrange multipliers:
    for ( Mecable * mec : objs )
    {
        real * f = vFOR + DIM * mec->matIndex();
        mec->getForces(f);
        mec->computeTensions(f);
    }
}


/**
This preforms:
 
    vFOR <- vFOR + Noise
    vRHS <- time_step * P * vFOR:
 
 Also prepare Projection diff is requested

 'rhs' and 'fff' are output. Input 'rnd' is a set of independent random numbers
*/
real browian1(Mecable* mec, real const* rnd, real alpha, real* fff, real time_step, real* rhs)
{
    real n = mec->addBrownianForces(rnd, alpha, fff);
    
    // Calculate the right-hand-side of the system in vRHS:
    mec->setSpeedsFromForces(fff, time_step, rhs);
    
    /*
     In this case, `fff` contains the true force in each vertex of the system
     and the Lagrange multipliers will represent the tension in the filaments
     */
    mec->storeTensions(fff);
    
#if ADD_PROJECTION_DIFF
    mec->makeProjectionDiff(fff);
#endif
    
    return n;
}


/**
 This solves the equation:
 
     ( Xnew - Xold ) / time_step = P * Force + Noise
 
 Explicit integration would lead to:
 
     Xnew = Xold + time_step * P * Force + Noise
 
 For implicit integration, we use a linearization of the force:
 
     Force(X) = M * X + B
 
 where M is a matrix and B is a vector, leading to:
 
     ( I - time_step * P * M ) ( Xnew - Xold ) = time_step * P * F + Noise
 
 where:
 
     F = M * Xold + B
 
     Noise = sqrt(2*kT*time_step*mobility) * Gaussian(0,1)
 
 Implicit integration ensures numerical stability. In the code,
 the matrix M is decomposed as:
 
     M = mB + mC + the rigidity terms of Mecables
 
 and the vectors are:
 
     'vPTS' is Xold
     'vBAS' is B
     'vRND' is Noise, a vector of calibrated Brownian terms
     'vFOR' is the force ( M * Xold + B ) and then ( M * Xnew + B )
     'vRHS' is the right-hand-side of the system (time_step * P * F + vRND)
     'vSOL' is `Xnew - Xold`, the solution to the system
 
 */
void Meca::solve(SimulProp const* prop, const int precond)
{
    assert_true( time_step == prop->time_step );
    
    prepareMatrices();
    
    // calculate forces before constraints in vFOR:
    calculateForces(vPTS, vBAS, vFOR);
    
#if ( DIM > 1 )
    addAllRigidity(vPTS, vFOR);
#endif
 
    /* 
     Fill `vRND` with Gaussian random numbers 
     This operation can be done in parallel, in a separate thread
     */
    RNG.gauss_set(vRND, dimension());
    
    /*
     Add Brownian motions to 'vFOR', and calculate vRHS by multiplying by mobilities.
     As Brownian terms are added, we record the magnitude of the typical smallest
     scalar contribution in `noiseLevel`. The dynamics will later be solved with 
     a residual that is proportional to this level:
     SimulProp::tolerance * noiseLevel
     As long as SimulProp::tolerance is smaller than 1, this should allow for a
     level of numerical error is small with respect to the Brownian noise in
     the system, and the results should be physically appropriate.
     */
    
    real noiseLevel = INFINITY;
    
    /*
     Add Brownian contributions and calculate Minimum value of it
      vFOR <- vFOR + Noise
      vRHS <- time_step * P * vFOR:
     */
#if NUM_THREADS > 1
    #pragma omp parallel num_threads(NUM_THREADS)
    {
        real local_res = INFINITY;
        Mecable ** mci = objs.begin() + omp_get_thread_num();
        while ( mci < objs.end() )
        {
            const index_t inx = DIM * (*mci)->matIndex();
            real n = browian1(*mci, vRND+inx, prop->kT/time_step, vFOR+inx, time_step, vRHS+inx);
            local_res = std::min(local_res, n);
            mci += NUM_THREADS;
        }
        //printf("thread %i min: %f\n", omp_get_thread_num(), local_res);
    #pragma omp critical
        noiseLevel = std::min(noiseLevel, local_res);
    }
#else
    for ( Mecable * mec : objs )
    {
        const index_t inx = DIM * mec->matIndex();
        real n = browian1(mec, vRND+inx, prop->kT/time_step, vFOR+inx, time_step, vRHS+inx);
        noiseLevel = std::min(noiseLevel, n);
    }
#endif
    
    if ( noiseLevel > 0 )
        noiseLevel *= time_step;
    else
        noiseLevel = 1.0;
    
    //printf("noiseLeveld = %8.2e   variance(vRHS) / estimate = %8.4f\n",
    //       noiseLevel, blas::nrm2(dimension(), vRHS) / (noiseLevel * sqrt(dimension())) );

#if NEW_CYTOPLASMIC_FLOW
    /**
     Includes a constant fluid flow displacing all the objects along
     */
    if ( prop->flow.norm() > REAL_EPSILON )
    {
        PRINT_ONCE("NEW_CYTOPLASMIC_FLOW code enabled\n");
        Vector flow_dt = prop->flow * time_step;
        for ( int p = 0; p < dimension(); ++p )
            flow_dt.add_to(vRHS+DIM*p);
    }
#endif
    
#if EXPLICIT_INTEGRATION
    /*
     This implements the forward Euler integration, for testing purposes.
     The result is very inefficient, since we have built the entire stiffness matrix,
     which is not necessary for this simple explicit scheme.
     */
    blas::xaxpy(dimension(), 1.0, vRHS, 1, vPTS, 1);
    
    for ( Mecable * mec : objs )
    {
        mec->getPoints(vPTS+DIM*mec->matIndex());
        mec->getForces(vFOR+DIM*mec->matIndex());
    }
    return;
#endif
    
    /*
     Choose the initial guess for the solution of the system (Xnew - Xold):
     we could use the solution at the previous step, or a vector of zeros.
     Using the previous solution could be advantageous if the speed were 
     somehow continuous. However, the system is without inertia. In addition,
     objects are considered in a random order to build the linear system, such
     that the blocks from two consecutive iterations do not match.
     From this, using zero for the initial guess seems safer:
     */
    zero_real(dimension(), vSOL);

    /*
     We now solve the system MAT * vSOL = vRHS  by an iterative method:
     the tolerance is in scaled to the contribution of Brownian
     motions contained in vRHS, assuming that the error is equally spread 
     along all degrees of freedom, this should work for tolerance << 1
     here a printf() can be used to check that the estimate is correct:
    */ 
    
    // NOTE: the tolerance to solve the system should be such that the solution
    // found does not depend on the initial guess.
    
    real residual = noiseLevel * prop->tolerance;

    /*
     With exact arithmetic, biConjugate Gradient should converge at most
     in a number of iterations equal to the size of the linear system, with
     each iteration involving 2 matrix-vector multiplications.
     We set here the max limit to the number of matrix-vector multiplication:
     */
    LinearSolvers::Monitor monitor(2*dimension(), residual);

    //------- call the iterative solver:

    if ( precond == 1 )
        computePreconditionner();

    //fprintf(stderr, "Solve precond %i size %6i\n", precond, dimension());

    // GMRES generally performs best:
    LinearSolvers::GMRES(*this, vRHS, vSOL, 64, monitor, allocator, mH, mV, temporary);

    //std::clog << "Solve size " << dimension() << "  precondition " << precond << "  " << residual << "\n";
    //fprintf(stderr, "    GMRES     count %4i  residual %10.6f\n", monitor.count(), monitor.residual());
#if ( 0 )
    // enable this to compare with another GMRES
    monitor.reset();
    zero_real(dimension(), vSOL);
    LinearSolvers::GMRES(*this, vRHS, vSOL, 32, monitor, allocator, mH, mV, temporary);
    fprintf(stderr, "    GMRES-32  count %4i  residual %10.6f\n", monitor.count(), monitor.residual());
#endif
#if ( 0 )
    // enable this to compare BCGS and GMRES
    monitor.reset();
    zero_real(dimension(), vSOL);
    if ( precond )
        LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
    else
        LinearSolvers::BCGS(*this, vRHS, vSOL, monitor, allocator);
    fprintf(stderr, "    BCGS     count %4i  residual %10.6f\n", monitor.count(), monitor.residual());
#endif
#if ( 0 )
    // enable this to compare with another implementation of biconjugate gradient stabilized
    monitor.reset();
    zero_real(dimension(), vSOL);
    LinearSolvers::bicgstab(*this, vRHS, vSOL, monitor, allocator);
    fprintf(stderr, "    bcgs  count %4i residual %10.6f\n", monitor.count(), monitor.residual());
#endif
    
    //------- in case the solver did not converge, we try other methods:
    
    if ( !monitor.converged() )
    {
        Cytosim::out("Solver failed: precond %i flag %i count %4i residual %.2e\n",
            precond, monitor.flag(), monitor.count(), monitor.residual());
        
        // try with different initial seed: vRHS
        monitor.reset();
        copy_real(dimension(), vRHS, vSOL);
        LinearSolvers::GMRES(*this, vRHS, vSOL, 255, monitor, allocator, mH, mV, temporary);
        Cytosim::out("     seed: count %4i residual %.2e\n", monitor.count(), monitor.residual());
        
        if ( !monitor.converged() )
        {
            monitor.reset();
            zero_real(dimension(), vSOL);
            
            if ( precond )
            {
                // try with the other method:
                LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                Cytosim::out("    BCGSP: count %4i residual %.2e\n", monitor.count(), monitor.residual());
            }
            else
            {
                // try with a preconditioner
                computePreconditionner();
                LinearSolvers::GMRES(*this, vRHS, vSOL, 127, monitor, allocator, mH, mV, temporary);
                Cytosim::out("    GMRES: count %4i residual %.2e\n", monitor.count(), monitor.residual());
                if ( !monitor.converged() )
                {
                    // try with other method:
                    monitor.reset();
                    zero_real(dimension(), vSOL);
                    LinearSolvers::BCGSP(*this, vRHS, vSOL, monitor, allocator);
                }
            }
            
            if ( !monitor.converged() )
            {
                // try with different GMRES parameters:
                monitor.reset();
                zero_real(dimension(), vSOL);
                LinearSolvers::GMRES(*this, vRHS, vSOL, 255, monitor, allocator, mH, mV, temporary);
                Cytosim::out("    GMRES(256): count %4i residual %.2e\n", monitor.count(), monitor.residual());
                
                if ( !monitor.converged() )
                {
                    // no method could converge... this is really bad!
                    std::stringstream oss;
                    oss << "Solve() failed to converge (" << monitor.count() << " iterations, " << monitor.residual() << ")";
                    throw Exception(oss.str());
                }
            }
        }
    }
    
#ifndef NDEBUG
    
    //check validity of the data:
    for( index_t ii = 0; ii < dimension(); ++ii )
    {
        if ( not_a_number(vSOL[ii]) )
        {
            fprintf(stderr, "Meca::solve produced NaN\n");
            abort();
        }
    }
    
#endif
    
    //add the solution of the system (=dPTS) to the points coordinates
    //blas::xaxpy(dimension(), 1.0, vSOL, 1, vPTS, 1);
    blas::add(dimension(), vSOL, vPTS);
    /*
     Re-calculate forces with the new coordinates, excluding bending elasticity
     and Brownian terms on the vertices.
     In this way the result returned to the fibers does not sum-up to zero,
     and is appropriate for example to calculate the effect of force on assembly.
     */
    calculateForces(vPTS, vBAS, vFOR);
    ready_ = 1;
    
    // report on the matrix type and size, sparsity, and the number of iterations
    if ( prop->verbose )
    {
        std::stringstream oss;
        oss << "Meca " << DIM << "x" << nbPts;
        oss << " block " << largestMecable();
        oss << " " << mB.what();
        if ( useMatrixC ) oss << " " << mC.what();
        oss << " precond " << precond;
        oss << " count " << monitor.count();
        //oss << " flag " << monitor.flag();
        oss << " residual " << monitor.residual() << "\n";
        if ( prop->verbose > 1 )
            std::clog << oss.str();
        else
            Cytosim::out << oss.str();
    }
}


// transfer newly calculated point coordinates back to Mecables
void Meca::apply()
{
    if ( ready_ )
    {
        ready_ = 0;
        
#if NUM_THREADS > 1
#pragma omp parallel num_threads(NUM_THREADS)
        {
            Mecable ** mci = objs.begin() + omp_get_thread_num();
            while ( mci < objs.end() )
            {
                Mecable * mec = *mci;
                mec->getForces(vFOR+DIM*mec->matIndex());
                mec->getPoints(vPTS+DIM*mec->matIndex());
                mci += NUM_THREADS;
            }
        }
#else
        for ( Mecable * mec : objs )
        {
            mec->getForces(vFOR+DIM*mec->matIndex());
            mec->getPoints(vPTS+DIM*mec->matIndex());
        }
#endif
    }
    else
    {
        //write(2, "extra Meca::apply() calls\n", 26);
    }
}


//------------------------------------------------------------------------------
#pragma mark - Debug/Output Functions

/**
 Extract the full matrix associated with matVect, in `mat` of size `sz`.
 */
void Meca::getSystem(index_t dim, real * mat) const
{
    if ( dim != dimension() )
        throw InvalidIO("wrong matrix dimension");
    
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    zero_real(dim, res);
    
    for ( index_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        multiply(src, res);
        blas::xcopy(dim, res, 1, mat+ii*dim, 1);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the full matrix associated with multiply(), in binary format
 */
void Meca::dumpSystem(FILE * file) const
{
    const index_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( index_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        multiply(src, res);
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the elasticity matrix
 */
void Meca::dumpElasticity(FILE * file) const
{
    const index_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( index_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        
        calculateForces(src, nullptr, res);
        
#if ( DIM > 1 )
        addAllRigidity(src, res);
#endif

#if ADD_PROJECTION_DIFF
        for ( Mecable const* mec : objs )
        {
            const index_t inx = DIM * mec->matIndex();
            if ( mec->hasProjectionDiff() )
                mec->addProjectionDiff(src+inx, res+inx);
        }
#endif
        
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save the projection matrix multiplied by the mobility
 */
void Meca::dumpMobility(FILE * file) const
{
    const index_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( index_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        
        zero_real(dim, res);
        
        for ( Mecable const* mec : objs )
        {
            const index_t inx = DIM * mec->matIndex();
            // this includes the mobility, but not the time_step:
            mec->setSpeedsFromForces(src+inx, 1.0, res+inx);
        }
        
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


/**
 Save matrix associated with the preconditionner, in binary format
 */
void Meca::dumpPreconditionner(FILE * file) const
{
    const index_t dim = dimension();
    real * src = new_real(dim);
    real * res = new_real(dim);
    
    zero_real(dim, src);
    
    for ( index_t ii = 0; ii < dim; ++ii )
    {
        src[ii] = 1.0;
        precondition(src, res);
        fwrite(res, sizeof(real), dim, file);
        src[ii] = 0.0;
    }
    
    free_real(res);
    free_real(src);
}


void Meca::dumpObjectID(FILE * file) const
{
    real * vec = new_real(largestMecable());
    
    for ( size_t ii = 0; ii < objs.size(); ++ii )
    {
        const unsigned nbp = objs[ii]->nbPoints();
        for ( unsigned p=0; p < nbp; ++p )
            vec[p] = ii;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), nbp, file);
    }
    
    free_real(vec);
}


void Meca::dumpDrag(FILE * file) const
{
    real * vec = new_real(largestMecable());
    
    for ( Mecable const* mec : objs )
    {
        const unsigned nbp = mec->nbPoints();
        const real drag = mec->dragCoefficient() / nbp;
        for ( unsigned p=0; p < nbp; ++p )
            vec[p] = drag;
        for ( int d = 0; d < DIM; ++ d )
            fwrite(vec, sizeof(real), nbp, file);
    }
    
    free_real(vec);
}


/**
 This dump the total matrix and some vectors in binary files.
 
 This matlab code should read the output:
 
     ord = load('ord.txt');
     time_step = load('stp.txt');
     obj = fread(fopen('obj.bin'), ord, 'double');
     drg = fread(fopen('drg.bin'), ord, 'double');
     sys = fread(fopen('sys.bin'), [ord, ord], 'double');
     ela = fread(fopen('ela.bin'), [ord, ord], 'double');
     mob = fread(fopen('mob.bin'), [ord, ord], 'double');
     con = fread(fopen('con.bin'), [ord, ord], 'double');
     pts = fread(fopen('pts.bin'), ord, 'double');
     rhs = fread(fopen('rhs.bin'), ord, 'double');
     sol = fread(fopen('sol.bin'), ord, 'double');
 
 To display the matrices:

     imshow(abs(sys))
     imshow(abs(ela))
 
 You can then compare the results with matlab's own iterative method,
 and compare the result using a scatter plot:
 
     x = bicgstab(sys, rhs, 0.001, ord);
     plot(x, sol, '.');
 
 */
void Meca::dump() const
{
    FILE * f = fopen("ord.txt", "w");
    fprintf(f, "%u\n", dimension());
    fclose(f);
    
    f = fopen("stp.txt", "w");
    fprintf(f, "%f\n", time_step);
    fclose(f);
 
    f = fopen("drg.bin", "wb");
    dumpDrag(f);
    fclose(f);
    
    f = fopen("obj.bin", "wb");
    dumpObjectID(f);
    fclose(f);
    
    f = fopen("rhs.bin", "wb");
    fwrite(vRHS, sizeof(real), dimension(), f);
    fclose(f);
    
    f = fopen("sol.bin", "wb");
    fwrite(vSOL, sizeof(real), dimension(), f);
    fclose(f);
    
    f = fopen("pts.bin", "wb");
    fwrite(vPTS, sizeof(real), dimension(), f);
    fclose(f);
    
    f = fopen("sys.bin", "wb");
    dumpSystem(f);
    fclose(f);
    
    f = fopen("ela.bin", "wb");
    dumpElasticity(f);
    fclose(f);
    
    f = fopen("mob.bin", "wb");
    dumpMobility(f);
    fclose(f);
    
    f = fopen("con.bin", "wb");
    dumpPreconditionner(f);
    fclose(f);
}

