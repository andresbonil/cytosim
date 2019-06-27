// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MONITOR_H
#define MONITOR_H

#include "real.h"

/// Iterative Solvers
namespace LinearSolvers
{
    /// records the number of iterations, and the convergence
    class Monitor
    {
    private:
        
        int      flag_;
        
        unsigned iter_,  iterMax_, iterOld_;
        
        real     resid_, residMax_;
        
    public:
        
        /// set the maximum number of iterations, and the residual threshold
        Monitor(unsigned i, real r) { reset(); iterMax_ = i; residMax_ = r; }
        
        /// reset interation count and achieved residual
        void reset() { flag_ = 0; iter_ = 0; resid_ = INFINITY; iterOld_ = 0; }
        
        /// increment iteration count
        void operator ++() { ++iter_; }
        
        /// increment iteration count
        void operator +=(unsigned i) { iter_ += i; }
       
        /// value of return flag
        int flag()       const { return flag_; }
        
        /// set flag to `f`
        void flag(const int f) { flag_ = f; }

        /// iteration count
        unsigned count()  const { return iter_; }
        
        /// last achieved residual
        real residual()  const { return resid_; }
        
        /// true if achieve residual < residual threshold
        bool converged() const { return resid_ < residMax_; }
        
        /// check given residual and return true if threshold is achieved
        bool finished(real);
        
        /// calculate residual from `x` and return true if threshold is achieved
        bool finished(unsigned size, const real* x);
        
        /// calculate residual from `x` and set flag to `f`
        void finish(int f, unsigned size, const real* x) { flag_ = f; finished(size, x); }
    };
}

#endif

