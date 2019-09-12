// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matsym.h"
#include "assert_macro.h"
#include "cblas.h"


//------------------------------------------------------------------------------
MatrixSymmetric::MatrixSymmetric()
{
    allocated_ = 0;
    val        = nullptr;
    in_charge  = 1;
}

//------------------------------------------------------------------------------
void MatrixSymmetric::allocate(size_t alc)
{
    if ( alc > allocated_ )
    {
        allocated_ = alc;
        free_real(val);
        val = new_real(alc*alc);
    }
}

//------------------------------------------------------------------------------
void MatrixSymmetric::deallocate()
{
    if ( in_charge )
        free_real(val);
    allocated_ = 0;
    val = nullptr;
}

//------------------------------------------------------------------------------
void MatrixSymmetric::reset()
{
    for ( index_t i = 0; i < size_ * size_; ++i )
        val[i] = 0;
}

//------------------------------------------------------------------------------
void MatrixSymmetric::scale( real alpha )
{
    for ( index_t i = 0; i < size_ * size_; ++i )
        val[i] *= alpha;
}

//------------------------------------------------------------------------------
real& MatrixSymmetric::operator()( index_t x, index_t y)
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    return val[ std::max(x,y) + msLDD * std::min(x,y) ];
}

//------------------------------------------------------------------------------
real* MatrixSymmetric::addr( index_t x, index_t y) const
{
    assert_true( x < size_ );
    assert_true( y < size_ );
    return val + ( std::max(x,y) + msLDD * std::min(x,y) );
}

//------------------------------------------------------------------------------
bool MatrixSymmetric::nonZero() const
{
    return true;
}

//------------------------------------------------------------------------------
size_t MatrixSymmetric::nbElements(index_t start, index_t stop) const
{
    assert_true( start <= stop );
    assert_true( stop <= size_ );

    return size_ * ( stop - start );
}

//------------------------------------------------------------------------------
std::string MatrixSymmetric::what() const
{
    return "full-symmetric";
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAdd( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, val, size_, X, 1, 1.0, Y, 1);
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAddIso2D( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, val, size_, X+0, 2, 1.0, Y+0, 2);
    blas::xsymv('L', size_, 1.0, val, size_, X+1, 2, 1.0, Y+1, 2);
}

//------------------------------------------------------------------------------
void MatrixSymmetric::vecMulAddIso3D( const real* X, real* Y ) const
{
    blas::xsymv('L', size_, 1.0, val, size_, X+0, 3, 1.0, Y+0, 3);
    blas::xsymv('L', size_, 1.0, val, size_, X+1, 3, 1.0, Y+1, 3);
    blas::xsymv('L', size_, 1.0, val, size_, X+2, 3, 1.0, Y+2, 3);
}


