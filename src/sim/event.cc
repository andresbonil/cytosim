// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"
#include "parser.h"


void Event::reset(real time)
{
    nextEvent = time + RNG.exponential() / rate;
}


Event::Event()
: code(""), recurrent(false), rate(0), nextEvent(0)
{
}


Event::Event(real time, Glossary& opt)
{
    if (!opt.set(code, "code")) code = "";
    if (!opt.set(rate, "rate")) rate = 0;
    if (!opt.set(recurrent, "recurrent")) recurrent = 0;
    reset(time);
}


Event::~Event()
{
    //Cytosim::log("destroying Event %p\n", this);
}


/// stochastic firing at specified rate
void Event::step(Simul& sim)
{
    if ( recurrent || sim.time() > nextEvent )
    {
        sim.relax();
        do {
            nextEvent += RNG.exponential() / rate;
            Parser(sim, 1, 1, 1, 1, 1).evaluate(code, ", in event:code");
        } while ( sim.time() > nextEvent );
        sim.prepare();
    }
}


void Event::write(Outputter& out) const
{
}


void Event::read(Inputter& in, Simul& sim, ObjectTag tag)
{
}
