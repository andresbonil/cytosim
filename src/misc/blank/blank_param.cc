// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// Created by Francois Nedelec on 6/27/07.

#include "blank_param.h"
#include "glossary.h"
#include <cmath>
#include <iomanip>

void BlankParam::clear() 
{
    max       = 1000;
    diffusion = 1;
    time_step = 0.001;
    time      = 0;
    delay     = 32;
    repeat    = 1;
    config    = "config.cym";
}


void BlankParam::read(Glossary& glos) 
{
    glos.set(max,        "max");
    glos.set(diffusion,  "diffusion");
    glos.set(time_step,  "time_step");
    glos.set(time,       "time");
    glos.set(delay,      "delay");
    glos.set(repeat,     "repeat");
    glos.set(config,     "config");

    if ( time < 0 )
        throw InvalidParameter("time must be >= 0");
    
    if ( diffusion < 0 )
        throw InvalidParameter("diffusion must be >= 0");

    /**
     We want for one degree of freedom to fulfill `var(dx) = 2 D dt`
     And we use: dx = diffusion_dt * RNG.sreal()
     Since `sreal()` is uniformly distributed, its variance is 1/3,
     and we need `diffusion_dt^2 = 6 D dt`
     */
    diffusion_dt = sqrt( 6.0 * diffusion * time_step );
}

/// formatted output of one parameter
template<typename C>
static  void write_value(std::ostream& os, std::string const& name, C const& c)
{
    os << " " << std::left << std::setw(20) << name << " = " << c << ";" << std::endl;
}

void BlankParam::write(std::ostream& os) const
{
    write_value(os, "max",       max);
    write_value(os, "diffusion", diffusion);
    write_value(os, "time_step", time_step);
    write_value(os, "time",      time);
    write_value(os, "delay",     delay);
    write_value(os, "repeat",    repeat);
    write_value(os, "config",    config);
}

