// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "sim.h"
#include "exceptions.h"
#include "glossary.h"
#include "simul_prop.h"
#include "wanderer_prop.h"
#include "wanderer.h"
#include "simul.h"


Hand * WandererProp::newHand(HandMonitor* h) const
{
    return new Wanderer(this, h);
}


void WandererProp::clear()
{
    DigitProp::clear();
    diffusion = 0;

}


void WandererProp::read(Glossary& glos)
{
    DigitProp::read(glos);
    glos.set(diffusion, "diffusion");
}


void WandererProp::complete(Simul const& sim)
{
    DigitProp::complete(sim);
    diff_rate         = diffusion / (step_size*step_size);
    diff_rate_2       = diff_rate * 2;
    diff_rate_dt      = diff_rate * sim.prop->time_step;
    diff2_rate_dt     = diff_rate_dt*2;
    step_size_kT      = step_size/sim.prop->kT;
}


void WandererProp::completeStiffness(Simul const& sim, real stiffness)
{
    U_step_kT_2 = (stiffness*step_size*step_size) / (2*sim.prop->kT);
}


void WandererProp::write_values(std::ostream & os) const
{
    DigitProp::write_values(os);
    write_value(os, "diffusion",  diffusion);
}

