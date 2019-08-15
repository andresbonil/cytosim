// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "event_set.h"
#include "iowrapper.h"
#include "glossary.h"
#include "simul.h"


void EventSet::step()
{
    for ( Event * e=first(); e; e=e->next() )
        e->step(simul);
}


Property* EventSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    return nullptr;
}


Object * EventSet::newObject(const ObjectTag tag, unsigned num)
{
    Event * e = nullptr;
    if ( tag == Event::TAG )
        e = new Event();
    return e;
}


/**
 @defgroup NewEvent How to create an Event
 @ingroup NewObject

 Specify a new Event:
 
     new event NAME
     {
         code = CODE;
         rate = POSITIVE_REAL;
         recurrent = [0, 1];
     }
 */
ObjectList EventSet::newObjects(const std::string&, Glossary& opt)
{
    ObjectList res;
    Event * e = new Event(simul.time(), opt);
    res.push_back(e);
    return res;
}

