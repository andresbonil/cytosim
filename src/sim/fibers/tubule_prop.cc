// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "tubule_prop.h"
#include <cmath>
#include "sim.h"
#include "tubule.h"
#include "glossary.h"
#include "exceptions.h"
#include "property_list.h"
#include "simul_prop.h"


Fiber* TubuleProp::newFiber() const
{
    return new Tubule(this);
}


void TubuleProp::clear()
{
    FiberProp::clear();
    
    dynamic_model[0] = 0;
    dynamic_model[1] = 0;
    growing_force[0] = INFINITY;
    growing_force[1] = INFINITY;

    rebirth_rate[0]  = 0;
    rebirth_rate[1]  = 0;
    rebirth_prob[0]  = 0;
    rebirth_prob[1]  = 0;

    for ( unsigned int n=0; n<4; ++n )
    {
        dynamic_trans1[n] = 0;
        dynamic_trans2[n] = 0;
        dynamic_speed1[n] = 0;
        dynamic_speed2[n] = 0;
    }
}


void TubuleProp::read(Glossary& glos)
{
    FiberProp::read(glos);
    
    glos.set(dynamic_model, 2,  "dynamic_model");
    glos.set(growing_force, 2,  "growing_force");
    glos.set(rebirth_rate,  2,  "rebirth_rate");
    
    glos.set(dynamic_trans1, 4, "dynamic_trans1");
    glos.set(dynamic_trans2, 4, "dynamic_trans2");
    glos.set(dynamic_speed1, 4, "dynamic_speed1");
    glos.set(dynamic_speed2, 4, "dynamic_speed2");
    
#ifdef BACKWARD_COMPATIBILITY
    
    if ( glos.set(growing_force, 2, "dynamic_force") )
        Cytosim::warn << "fiber:dynamic_force was renamed growing_force" << std::endl;
    
    int f = 0;
    if ( glos.set(f, "fate", {{"none", 0}, {"destroy", 1}, {"rescue", 2}}) )
    {
        Cytosim::warn << "fiber:fate is deprecated: use `persistent` and `rebirth_rate`" << std::endl;
        persistent = ( f != 1 );
        rebirth_rate[0] = ( f == 2 ? INFINITY : 0 );
    }

#endif
}


void TubuleProp::complete(Simul const& sim)
{
    FiberProp::complete(sim);
    
    const real time_step = sim.prop->time_step;
    
    if ( growing_force[0] <= 0 )
        throw InvalidParameter("fiber:growing_force[0] should be > 0");
    if ( growing_force[1] <= 0 )
        throw InvalidParameter("fiber:growing_force[1] should be > 0");
    
    if ( rebirth_rate[0] < 0 )
        throw InvalidParameter("fiber:rebirth_rate[0] should be >= 0");
    if ( rebirth_rate[1] < 0 )
        throw InvalidParameter("fiber:rebirth_rate[1] should be >= 0");

    rebirth_prob[0] = 1 - exp( -rebirth_rate[0] * time_step );
    rebirth_prob[1] = 1 - exp( -rebirth_rate[1] * time_step );

    for ( unsigned int n=0; n<4; ++n )
    {
        if ( dynamic_trans1[n] < 0 )
            throw InvalidParameter("tubule:dynamic_trans1[",n,"] should be >= 0");
        if ( dynamic_trans1[n] * time_step > sim.prop->acceptable_rate )
            throw InvalidParameter("tubule:dynamic_trans1[",n,"] is too high: decrease time_step");
        
        if ( dynamic_trans2[n] < 0 )
            throw InvalidParameter("tubule:dynamic_trans2[",n,"] should be >= 0");
        if ( dynamic_trans2[n] * time_step > sim.prop->acceptable_rate )
            throw InvalidParameter("tubule:dynamic_trans2[",n,"] is too high: decrease time_step");
    }
}


void TubuleProp::write_values(std::ostream & os) const
{
    FiberProp::write_values(os);
    write_value(os, "dynamic_model",  dynamic_model, 2);
    write_value(os, "growing_force",  growing_force, 2);
    write_value(os, "rebirth_rate",   rebirth_rate, 2);
    write_value(os, "dynamic_speed1", dynamic_speed1, 4);
    write_value(os, "dynamic_speed2", dynamic_speed2, 4);
    write_value(os, "dynamic_trans1", dynamic_trans1, 4);
    write_value(os, "dynamic_trans2", dynamic_trans2, 4);
}

