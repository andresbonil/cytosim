// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "couple_long.h"
#include "couple_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

CoupleLong::CoupleLong(CoupleProp const* p, Vector const& w)
: Couple(p, w), mArm(nullTorque)
{
}


CoupleLong::~CoupleLong()
{
}

//------------------------------------------------------------------------------

Torque CoupleLong::calcArm(Interpolation const& pt, Vector const& pos, real len)
{
    Vector off = pt.pos1() - pos;
    if ( modulo )
        modulo->fold(off);
#if ( DIM >= 3 )
    off = cross(off, pt.diff());
    real n = off.norm();
    if ( n > REAL_EPSILON )
        return off * ( len / n );
    else
        return pt.diff().randOrthoU(len);
#else
    return std::copysign(len, cross(off, pt.diff()));
#endif
}

//------------------------------------------------------------------------------

/*
Note that, since `mArm` is calculated by setInteraction(),
the result of sidePos() will be incorrect if 'solve=0'
*/
Vector CoupleLong::sidePos() const
{
#if ( DIM == 1 )
    
    return cHand1->pos();
    
#elif ( DIM == 2 )
    
    return cHand1->pos() + cross(mArm, cHand1->dirFiber());
    
#else
    
    ///\todo: change formula to match interSideLink3D
    return cHand1->pos() + mArm;

#endif
}


/**
 Calculates the force for the interSideLink()
 */
Vector CoupleLong::force() const
{
    Vector d = cHand2->pos() - CoupleLong::sidePos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


/**
 This uses interSideLink().
 
 Another possibility would be SideSideLink, which is fully symmetric.
 */
void CoupleLong::setInteractions(Meca & meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
    //meca.addSideSideLink(pt1, pt2, prop->length, prop->stiffness);
    
    /*
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    
#if ( DIM == 2 )
    
    mArm = calcArm(pt1, pt2.pos(), prop->length);
    meca.addSideLink2D(pt1, pt2, mArm, prop->stiffness);
    
#elif ( DIM >= 3 )

    mArm = calcArm(pt1, pt2.pos(), prop->length);
    meca.addSideLink3D(pt1, pt2, mArm, prop->stiffness);
    //@todo Couple::setInteractions() use interSideLink3D()
    
#endif
    
}


