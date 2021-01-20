// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "crosslink_long.h"
#include "crosslink_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------

CrosslinkLong::CrosslinkLong(CrosslinkProp const* p, Vector const& w)
: Crosslink(p, w), mArm(nullTorque)
{
}


CrosslinkLong::~CrosslinkLong()
{
}

//------------------------------------------------------------------------------

Torque CrosslinkLong::calcArm(Interpolation const& pt, Vector const& pos, real len)
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
Vector CrosslinkLong::sidePos() const
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
Vector CrosslinkLong::force() const
{
    Vector d = cHand2->pos() - CrosslinkLong::sidePos();
    
    //correct for periodic space:
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


/**
 This uses interSideLink().
 
 Another possibility would be SideSideLink, which is fully symmetric.
 */
void CrosslinkLong::setInteractions(Meca & meca) const
{
    Interpolation const& pt1 = cHand1->interpolation();
    Interpolation const& pt2 = cHand2->interpolation();
    
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
    //@todo CrosslinkLong:setInteractions() use interSideLink3D()
    
#endif
    
    //meca.addSideSideLink( cHand1->interpolation(), cHand2->interpolation(), prop->length, prop->stiffness );
}


