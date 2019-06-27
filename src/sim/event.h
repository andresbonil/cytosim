// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef EVENT_H
#define EVENT_H

#include "assert_macro.h"
#include "object.h"

class Meca;
class Simul;
class Glossary;

/// Performs actions on the simulation
/** 
*/
class Event: public Object
{
    
    friend class EventSet;
    
public:
    
    /// code to be executed
    std::string code;
    
    /// true if event executes at every time step
    bool        recurrent;

    /// rate at which code is executed
    real        rate;
    
    /// time of next event
    real        nextEvent;
    
public:

    /// default constructor
    Event();
    
    /// constructor
    Event(real time, Glossary&);

    /// destructor
    virtual ~Event();
    
    /// initialize counters
    void reset(real time);
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 't';

    /// an ASCII character identifying the class of this object
    ObjectTag       tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return nullptr; }

    //--------------------------------------------------------------------------
    
    /// monte-carlo simulation step
    void      step(Simul&);
    
    /// add interactions to the Meca
    void      setInteractions(Meca &) const {}
    
    
    /// a static_cast<> of Node::next()
    Event *   next()  const  { return static_cast<Event*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    Event *   prev()  const  { return static_cast<Event*>(nPrev); }
    

    /// read
    void      read(Inputter&, Simul&, ObjectTag);
    
    /// write
    void      write(Outputter&) const;
};


#endif
