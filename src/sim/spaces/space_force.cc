// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_force.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


SpaceForce::SpaceForce(SpaceProp const* p)
: Space(p)
{
    force.reset();
    center.reset();
    stiffness = 0;
}

void SpaceForce::resize(Glossary& opt)
{
    opt.set(force, "force");
    opt.set(center, "center");
    opt.set(stiffness, "stiffness");
}


void SpaceForce::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-1, -1, -1);
    sup.set( 1,  1,  1);
}


real SpaceForce::volume() const
{
    throw InvalidParameter("invalid use of space `force'");
    return -1;
}


Vector SpaceForce::project(Vector const& w) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const
{
    throw InvalidParameter("Invalid use of space `force'");
}


void SpaceForce::setInteractions(Meca & meca, FiberSet const&) const
{
    real *const base = meca.base();

    if ( stiffness > 0 )
    {
        Vector sc = stiffness * center;
        const unsigned nbp = meca.nb_points();
        // generate an isotropic squeezing force:
        for ( unsigned p = 0; p < nbp; ++p )
        {
            meca.mB(p,p) -= stiffness;
            sc.add_to(base+DIM*p);
        }
    }
    else
    {
        const unsigned nbu = meca.dimension();
        for ( unsigned u = 0; u < nbu; u += DIM )
            force.add_to(base+u);
    }
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceForce::draw() const
{
    gle::gleArrow(center, center+force, 0.1 * force.norm());
    return true;
}

#else

bool SpaceForce::draw() const
{
    return false;
}


#endif


