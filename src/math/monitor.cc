// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "monitor.h"
#include <cstdio>
#include "cblas.h"

/**
 accept given residual and test against threshold
 */
bool LinearSolvers::Monitor::finished(real resid)
{
    resid_ = resid;

    if ( iter_ > iterMax_ )
        return true;

    return ( resid < residMax_ );
}


/**
 Here is defined which norm is used to measure the residual
 */
bool LinearSolvers::Monitor::finished(const unsigned size, const real * x)
{
    //fprintf(stderr, "Solver %3u residual %9.6f %9.6f\n", iter_, blas::nrm2(size, x), blas::nrm8(size, x));

#if ( 1 )
    // use the 'infinite' norm (i.e. the largest element)
    real resid = blas::nrm8(size, x);
#else
    // use the standard Euclidian norm:
    real resid = blas::nrm2(size, x);
#endif
#if ( 1 )
    if ( iter_ > iterOld_+63 )
    {
        if ( resid > 2*resid_ )
        {
            printf("Warning: slow convergence (time_step may be too big)");
            printf(" residual %9.6f at iteration %3u, %9.6f at %3u\n", resid_, iterOld_, resid, iter_);
        }
        iterOld_ = iter_;
    }
#endif
    
    if ( resid != resid )
    {
        fprintf(stderr, "Solver diverged at step %3u (residual is not-a-number)\n", iter_);
        return true;
    }
    
    return finished(resid);
}


