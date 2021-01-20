// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "wrist_long.h"
#include "simul.h"
#include "meca.h"
#include "modulo.h"


extern Modulo const* modulo;


WristLong::WristLong(SingleProp const* sp, Mecable const* mec, const unsigned pti)
: Wrist(sp, mec, pti), mArm(nullTorque)
{
}


WristLong::~WristLong()
{
}

//------------------------------------------------------------------------------

Torque WristLong::calcArm(Interpolation const& pt, Vector const& pos, real len)
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
Vector WristLong::force() const
{
    assert_true( sHand->attached() );
    Vector d = posFoot() - WristLong::sidePos();
    
    if ( modulo )
        modulo->fold(d);
    
    return prop->stiffness * d;
}


//------------------------------------------------------------------------------

/*
Note that, since `mArm` is calculated by setInteraction(),
the result of sidePos() will be incorrect if 'solve=0'
*/
Vector WristLong::sidePos() const
{
#if ( DIM > 1 )
    return sHand->pos() + cross(mArm, sHand->dirFiber());
#endif
    return sHand->pos();
}

/**
 Using a Meca::interSideLink()
 */
void WristLong::setInteractions(Meca & meca) const
{
    Interpolation const& pt = sHand->interpolation();
    
    /* 
     The 'arm' is recalculated each time, but in 2D at least,
     this maybe not necessary, as switching should be rare.
     */
    
#if ( DIM == 2 )
    
    mArm = calcArm(pt, posFoot(), prop->length);
    if ( anchor.rank() == 1 )
        meca.addSideLink2D(pt, anchor.vertex0(), mArm, prop->stiffness);
    else
        throw InvalidParameter("unfinished WristLong::setInteractions(length>0, Interpolation4)");

#elif ( DIM >= 3 )
    
    mArm = calcArm(pt, posFoot(), prop->length);
    if ( anchor.rank() == 1 )
        meca.addSideLink3D(pt, anchor.point(), mArm, prop->stiffness);
    else
        throw InvalidParameter("unfinished WristLong::setInteractions(length>0, Interpolation4)");
    
#endif
}

