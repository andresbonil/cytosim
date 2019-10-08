// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MONITOR_H
#define MONITOR_H

#include "real.h"
#include "cblas.h"

/// Iterative Solvers
namespace LinearSolvers
{
    /// records the number of iterations, and the convergence
    class Monitor
    {
    private:
        
        /// flag
        int      flag_;
        
        /// counter for iterations or number of matrix-vector operations
        unsigned cnt_,  cntMax_, cntOld_;
        
        /// desired residual
        real     resMax_;

        /// achieved residual
        real     res_;
        
    public:
        
        /// set the maximum number of iterations, and the residual threshold
        Monitor(unsigned i, real r) { reset(); cntMax_ = i; resMax_ = r; }
        
        /// reset state variables (counters, flags and residual)
        void reset() { flag_ = 0; cnt_ = 0; res_ = INFINITY; cntOld_ = 32; }
        
        /// increment counter
        void operator ++() { ++cnt_; }
        
        /// increment counter by `i`
        void operator +=(unsigned i) { cnt_ += i; }
       
        /// value of return flag
        int flag()       const { return flag_; }
        
        /// set flag to `f`
        void flag(const int f) { flag_ = f; }

        /// iteration count
        unsigned count() const { return cnt_; }
        
        /// last achieved residual
        real residual()  const { return res_; }
        
        /// true if achieve residual < residual threshold
        bool converged() const { return res_ < resMax_; }
        
        /// check given residual and return true if threshold is achieved
        bool finished(real res)
        {
            res_ = res;
            
            if ( cnt_ > cntMax_ )
                return true;
            
            return ( res < resMax_ );
        }

        /// calculate residual from `x` and return true if threshold is achieved
        bool finished(unsigned size, const real* x)
        {
            //fprintf(stderr, "Solver %3u residual %9.6f %9.6f\n", cnt_, blas::nrm2(size, x), blas::nrm8(size, x));
            
#if ( 1 )
            // use the 'infinite' norm (i.e. the largest element)
            real res = blas::nrm8(size, x);
#else
            // use the standard Euclidian norm:
            real res = blas::nrm2(size, x);
#endif
#if ( 1 )
            if ( cnt_ > cntOld_+128 )
            {
                if ( res > 2*res_ )
                {
                    printf("Warning: slow convergence (time_step may be too big)");
                    printf(" residual %.3e at iteration %3u, %.3e at %3u\n", res_, cntOld_, res, cnt_);
                }
                cntOld_ = cnt_;
            }
#endif
            
            if ( res != res )
            {
                fprintf(stderr, "Solver diverged at step %3u (residual is not a number)\n", cnt_);
                return true;
            }
            
            return finished(res);
        }

        /// calculate residual from `x` and set flag to `f`
        void finish(int f, unsigned size, const real* x) { flag_ = f; finished(size, x); }
    };
}

#endif

