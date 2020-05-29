// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <cmath>
#include "growing_fiber_prop.h"
#include "growing_fiber.h"
#include "property_list.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul.h"


Fiber* GrowingFiberProp::newFiber() const
{
    return new GrowingFiber(this);
}


void GrowingFiberProp::clear()
{
    FiberProp::clear();
    
    for ( int i = 0; i < 2; ++i )
    {
        growing_speed[i]     = 0;
        growing_off_speed[i] = 0;
        growing_force[i]     = INFINITY;
        shrink_outside[i]    = false;
        shrinking_speed[i]   = 0;
        shrinking_space = "";
    }
}


void GrowingFiberProp::read(Glossary& glos)
{
    FiberProp::read(glos);
    
    glos.set(growing_force,     2, "growing_force");
    glos.set(growing_speed,     2, "growing_speed");
    glos.set(growing_off_speed, 2, "growing_off_speed");
    glos.set(shrink_outside,    2, "shrink_outside");
    glos.set(shrinking_speed,   2, "shrinking_speed");
    glos.set(shrinking_space, "shrinking_space");
}


void GrowingFiberProp::complete(Simul const& sim)
{
    FiberProp::complete(sim);

    for ( int i = 0; i < 2; ++i )
    {
        if ( growing_force[i] <= 0 )
            throw InvalidParameter("fiber:growing_force should be > 0");
        
        if ( growing_speed[i] < 0 )
            throw InvalidParameter("fiber:growing_speed should be >= 0");
        
        if ( shrinking_speed[i] > 0 )
            throw InvalidParameter("fiber:shrinking_speed should be <= 0");

        growing_speed_dt[i] = growing_speed[i] * sim.time_step();
        growing_off_speed_dt[i] = growing_off_speed[i] * sim.time_step();
        shrinking_speed_dt[i] = shrinking_speed[i] * sim.time_step();
    }
    if (shrink_outside[0]||shrink_outside[1])
    {
        if (shrinking_space!="")
            shrinking_space_ptr = sim.findSpace(shrinking_space);
        else
        {
            shrinking_space_ptr = confine_space_ptr;
            shrinking_space = confine_space;
        }
    }
}


void GrowingFiberProp::write_values(std::ostream& os) const
{
    FiberProp::write_values(os);
    write_value(os, "growing_speed",     growing_speed, 2);
    write_value(os, "growing_off_speed", growing_off_speed, 2);
    write_value(os, "growing_force",     growing_force, 2);
    write_value(os, "shrink_outside",    shrink_outside, 2);
    write_value(os, "shrinking_speed",   shrinking_speed, 2);
    if (shrink_outside[0]||shrink_outside[1])
        write_value(os, "shrinking_space",   shrinking_space);
    
}

