// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef MECAFIL_H
#define MECAFIL_H

#include "chain.h"
#include "fiber_prop.h"  // needed for NEW_FIBER_LOOP

class Matrix;

/// incompressible Chain with bending elasticity
/**
 Implements the methods of a Mecable for the Chain:
 
 -# setSpeedsFromForces() includes longitudinal incompressibility,
 which means keeping successive points equidistants:
 norm( point(p+1) - point(p) ) = segmentation()
 
 -# addRigidity() implements bending elasticity.
 .
*/
class Mecafil : public Chain
{
private:
    
    /// Lagrange multipliers associated with longitudinal imcompressibility
    real   *    rfLag;
    
    /// stored normalized differences of successive vertices (array of size DIM*nbSegments)
    real   *    rfDiff;
    
    /// memory allocated to hold nbPoints() values (used as temporary variables)
    real   *    rfLLG, * rfVTP;
    
    /// J*J', a nbSegments^2 matrix. We store the diagonal and one off-diagonal
    real   *    mtJJt, * mtJJtU;
    
    /// vector for the projection correction of size nbSegments
    real   *    mtJJtiJforce;
    
    /// true if all elements of mtJJtiJforce[] are null
    bool        useProjectionDiff;
    
protected:
    
    /// mobility of the points (all points have the same drag coefficient)
    real        rfDragPoint;
    
    /// rigidity scaling factor used in addRigidity()
    real        rfRigidity;
    
#if NEW_FIBER_LOOP
    /// link filament into a loop
    bool        rfRigidityLoop;
#endif
    
    /// calculate the normalized difference of successive vertices in rfDiff[]
    void        storeDirections();

private:
    
    /// reset the memory pointers for the projection
    void        buildProjection();
    
    /// allocate memory for the projection
    void        allocateProjection(size_t);
    
    /// free the memory for the projection
    void        destroyProjection();

public:
    
    /// Constructor
    Mecafil();
    
    /// copy constructor
    Mecafil(Mecafil const&);
    
    /// copy assignment
    Mecafil& operator=(Mecafil const&);
    
    /// Destructor
    virtual    ~Mecafil();
    
    /// allocate memory
    size_t      allocateMecable(size_t);
    
    /// free allocated memory
    void        release();

    /// compute Lagrange multiplier corresponding to the longitudinal tensions in the segments
    void        computeTensions(const real* force);
    
    /// save Lagrange multipliers computed in setSpeedsFromForces()
    void        storeTensions(const real* force);

    /// longitudinal force along segment `p`
    /**
     Tensions are calculated as the Lagrange multipliers associated with the
     constrains of conserved segments lengths.
     The tension is:
     - positive when the fiber is being pulled
     - negative when the fiber is under compression
     .
     */
    real        tension(unsigned p) const { assert_true(p+1<nPoints); return rfLag[p]; }
    
    /// total drag-coefficient of object (force = drag * speed)
    real        dragCoefficient() const { return  nPoints * rfDragPoint; }
    
    /// drag coefficient of one point
    real        leftoverDrag() const { return rfDragPoint; }

    //--------------------- Projection  / Dynamics
    
    /// prepare for projection
    void        makeProjection();
    
    /// prepare the correction to the projection
    void        makeProjectionDiff(const real* );
    
    /// add the contribution from the projection correction
    void        addProjectionDiff(const real*, real*) const;
    
    /// true if addProjectionDiff() does something
    bool        hasProjectionDiff() const { return useProjectionDiff; }

    /// add displacements due to the Brownian motion to rhs[]
    real        addBrownianForces(real const* rnd, real alpha, real* rhs) const;
    
    /// calculate the speeds from the forces, including projection
    void        setSpeedsFromForces(const real* X, real alpha, real* Y) const;
    
    /// print projection matrix
    void        printProjection(std::ostream&) const;

    //--------------------- Rigidity

    /// add the rigidity force corresponding to configuration X into vector Y
    void        addRigidity(const real* X, real* Y) const;
    
    /// add rigidity terms to upper side of matrix
    void        addRigidityUpper(real*, unsigned) const;
    
    /// add rigidity terms on three specified points
    void        addRigidity(const real* X, real* Y, unsigned, unsigned, unsigned) const;

};


#endif
