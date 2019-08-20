// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "monitor.h"
#include <cstdio>
#include "cblas.h"

/**
 accept given residual and test against threshold
 */
bool LinearSolvers::Monitor::finished(real resid)
{
#if ( 0 )
    if ( iter_ > 1024 && iter_ % 16 == 0 )
        printf("Monitor:  %4i residual %12.6f\n", iter_, resid);
#endif
    
#if ( 0 )
    if ( iter_ > iterOld_+63 )
    {
        if ( resid >= resid_ )
        {
            printf("Warning: slow convergence (time_step may be too big)");
            printf(" residual %9.6f at iteration %3i, %9.6f at %3i\n", resid_, iterOld_, resid, iter_);
        }
        iterOld_ = iter_;
    }
#endif
    
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
    // use the 'infinite' norm
    //real resid = blas::nrm8(size, x);
    
    // use the standard Euclidian norm:
    real resid = blas::nrm2(size, x);

#if ( 1 )
    if ( iter_ > iterOld_+63 )
    {
        if ( resid >= 2*resid_ )
        {
            printf("Warning: slow convergence (time_step may be too big)");
            printf(" residual %9.6f at iteration %3u, %9.6f at %3u\n", resid_, iterOld_, resid, iter_);
        }
        iterOld_ = iter_;
    }
#endif
    
#if ( 0 )
    std::cerr << "Solver step " << iter_ << " residual " << resid << '\n';
#endif
    
    if ( resid != resid )
    {
        fprintf(stderr, "Solver diverged at step %3u residual is not-a-number\n", iter_);
        return true;
    }
    
    return finished(resid);
}


