// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef LINEAR_OPERATOR_H
#define LINEAR_OPERATOR_H

#include "real.h"
#include "cblas.h"

namespace LinearSolvers
{
    /// defines the functions that defines the linear transformation
    class LinearOperator
    {
    public:
        /// size of the matrix M
        virtual unsigned dimension() const = 0;
        
        /// apply operator to a vector ( Y <- M * X )
        virtual void multiply(const real* X, real* Y) const = 0;
        
        /// apply transposed operator to vector ( Y <- transpose(M) * X )
        virtual void trans_multiply(const real* X, real* Y) const {}
        
        /// apply preconditionning ( Y <- P * X )
        virtual void precondition(const real* X, real* Y) const
        {
            copy_real(dimension(), X, Y);
        }
    };
}

#endif

