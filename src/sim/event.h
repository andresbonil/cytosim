// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef EVENT_H
#define EVENT_H

#include "object.h"

class Meca;
class Simul;
class Glossary;

/// an Event performs action on the simulation by executing code
/**
 An Event is a class that can perform some action in the simulation world,
 at regular interval or at stochastic time with a specified rate.
 The action is specified code interpreted by cytosim's parser.
 This can be used for example to add or remove objects.
 
 It is a special class that is not associated with a Property,
 and can be created with 'new' without a preceeding 'set'.
 
 Events are not saved to trajectory files.
*/
class Event: public Object
{
    
    friend class EventSet;
    
    /// clear member variables
    void clear();
    
public:
    
    /**
     @defgroup EventPar Parameters of Event
     @ingroup Parameters
     These are the parameters for an Event
     @{
     */

    /// code to be executed
    std::string activity;
    
    /// true if event executes at every time step
    bool        recurrent;

    /// rate at which code is executed
    real        rate;
    
    ///@}
    
    /// time of next event
    real        nextTime;
    
public:

    /// default constructor
    Event() { clear(); }
    
    /// constructor
    Event(real time, Glossary&);

    /// destructor
    virtual ~Event();
    
    /// initialize counters
    void reset(real time);
    
    /// a unique character identifying the class
    static const ObjectTag TAG = 'q';

    /// an ASCII character identifying the class of this object
    ObjectTag       tag() const { return TAG; }

    /// returns 0, since Event have no Property
    Property const* property() const { return nullptr; }

    //--------------------------------------------------------------------------
    
    /// monte-carlo simulation step
    void      step(Simul&);
    
    /// add interactions to a Meca
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
