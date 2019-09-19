// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "crosslink.h"
#include "crosslink_prop.h"
#include "exceptions.h"
#include "random.h"
#include "modulo.h"
#include "space.h"
#include "meca.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
Crosslink::Crosslink(CrosslinkProp const* p, Vector const& w)
: Couple(p, w), prop(p)
{
}


Crosslink::~Crosslink()
{
    prop = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Simulates:
 - diffusive motion
 - attachment
 .
 */
void Crosslink::stepFF(const FiberGrid& grid)
{
    diffuse();
    
    // confinement:
    if ( !prop->confine_space_ptr->inside(cPos) )
        cPos = prop->confine_space_ptr->bounce(cPos);
    
    if ( modulo )
        modulo->fold(cPos);
    
    // activity:
    cHand1->stepUnattached(grid, cPos);
    cHand2->stepUnattached(grid, cPos);
}


void Crosslink::setInteractions(Meca & meca) const
{
    assert_true( attached1() && attached2() );
    
    meca.addLink(cHand1->interpolation(), cHand2->interpolation(), prop->stiffness);
}

