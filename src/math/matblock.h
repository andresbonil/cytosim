// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MATBLOCK_H
#define MATBLOCK_H

#include "real.h"
#include "matrix.h"
#include "assert_macro.h"
#include <cstdio>
#include <string>

/// A non-symmetric real Matrix, diagonal by blocks.
/**
 (not used in Cytosim)
 */
class MatrixOfBlocks
{    
private:
    
    /// size of matrix
    index_t    size_;

    /// size of memory which has been allocated
    size_t     allocated_;
    
    /// number of blocks
    index_t    block_cnt;
    
    /// array of pointers to the blocks
    real**     block_;
    
    /// array containing the size of each block
    index_t*   block_size;
    
    /// array specifying the allocated size of each block
    size_t*    block_alc;
        
public:
    
    /// return the size of the matrix
    index_t size() const { return size_; }
    
    /// change the size of the matrix
    void resize(index_t s) { allocate(s); size_=s; }
    
    /// default constructor
    MatrixOfBlocks();
    
    /// the deallocation
    void deallocate();
    
    /// allocate the matrix to be able to hold nb_block (arg 1) blocks
    void allocate(size_t nb_block);
    
    /// allocate block b (arg 1) to be capable of holding (size*size) (arg 2)
    void setBlockSize( unsigned int b, unsigned int size);
    
    /// default destructor
    virtual ~MatrixOfBlocks() { deallocate(); }
    
    
    /// return the address of first element in block ii
    real* block( const unsigned int ii ) const
    {
        assert_true( block_ );
        assert_true( ii < block_cnt );
        assert_true( block_[ii] );
        assert_true( block_size[ii] <= block_alc[ii] );
        return block_[ii];
    }
    
    /// returns the number of blocks
    unsigned int nbBlocks() const 
    {
        return block_cnt;
    }
    
    /// calculate and return the total size of the matrix
    unsigned int calculateSize();
        
    /// returns the size of block ii
    unsigned int blockSize( const unsigned int ii ) const 
    {
        assert_true( block_ );
        assert_true( ii < block_cnt );
        return block_size[ii];
    }
    
    /// the address holding element (ii, jj)
    real* addr( index_t ii, index_t jj ) const;
    
    /// returns the address of element at (x, y), allocating if necessary
    virtual real& operator()( index_t x, index_t y );
    
    /// reset all the values in block ii
    void setBlockToZero( unsigned int ii );
    
    /// reset the entire matrix
    void reset();
    
    /// scale block ii
    void scaleBlock( unsigned int ii, real a );
    
    /// scale the entire matrix
    void scale( real a );
    
    /// total number of elements allocated
    size_t nbElements() const;
    
    /// size of the biggest block
    size_t maxBlockSize() const;
    
    /// vector multiplication: Y <- M * X
    void vecMulAdd( const real* X, real* Y) const;
    
    /// vector multiplication: Y <- M * X
    void vecMul( const real* X, real* Y) const;
    
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 2 * size(M)
    void vecMulAddIso2D( const real* X, real* Y ) const { }
    
    /// isotropic vector multiplication: Y = Y + M * X, size(X) = size(Y) = 3 * size(M)
    void vecMulAddIso3D( const real* X, real* Y ) const { }
    
    /// maximum of the absolute value of all elements
    real norm_inf() const;
    
    /// printf debug function in sparse mode: i, j : value
    void printSparse(std::ostream&) const;
    
    /// returns a string which a description of the type of matrix
    std::string what() const;
};


#endif
