// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef MODULO_H
#define MODULO_H

#include "real.h"
#include "vector.h"

/// Modulo implements periodic boundary conditions
/**
 This class is used to apply periodic boundaries conditions in one or 
 in multiple dimensions in space, X Y or Z.

 The methods are not virtual to avoid the calling overload in C++.
 */
class Modulo
{
private:
    
    /// bitfield indicating the dimensions that are periodic
    int   mMode;

    /// half-period in each dimension
    real  mSize[DIM];
    
    /// adjust 'x' to canonical image within periodicity 'p':
    static inline void fold(real& x, const real p)
    {
        if ( fabs(x) > p )
        {
            real i = std::copysign(p, x);
            do
                x = x - 2.0 * i;
            while ( fabs(x) > p );
        }
    }

public:
    
    /// constructor
    Modulo() { mMode = 0; for (int d=0; d<DIM; ++d) mSize[d] = 0; }

    /// destructor
    ~Modulo() {}
    
    /// disable periodicity in all dimensions
    void disable() { mMode = 0; }
    
    /// enable periodicity in dimension 'd'
    void enable(int d, real size);
    
    /// true if at least one direction has periodic boundaries
    bool isPeriodic() const { return mMode; }

    /// true if direction `d` has periodic boundaries
    bool isPeriodic(int d) const { return mMode & (1<<d); }
    
    /// return the d-th direction of periodicity
    const Vector periodicity(int d) const;
    
    /// shift `pos` to its canonical image, which is the one closest to the origin
    void         fold(Vector & pos) const;
    
    /// shift `pos` to its image which is closest to `ref`
    void         fold(Vector & pos, Vector const& ref) const;
    
    /// return translation necessary to bring `pos` to its canonical image
    const Vector offset(Vector const& pos) const;
    
    /// set `pos` to its canonical image, and return offset = pos - fold(pos)
    void         foldOffset(Vector& pos, Vector& off) const;
    
};

#endif


