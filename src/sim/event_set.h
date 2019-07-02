// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef EVENT_SET_H
#define EVENT_SET_H

#include "object_set.h"
#include "event.h"
class Simul;


/// a list of Event
/**
 */
class EventSet : public ObjectSet
{
public:
    
    /// creator
    EventSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "event"; }
    
    /// create a new property of category `cat` for a class `name`
    Property *  newProperty(const std::string& cat, const std::string& name, Glossary&) const;
    
    /// create objects of class `name`, given the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(ObjectTag, unsigned);
    
    //--------------------------
    
    /// first object
    Event *     first() const
    {
        return static_cast<Event*>(nodes.front());
    }
    
    ///  return pointer to the Object of given ID, or zero if not found
    Event *     findID(ObjectID n) const
    {
        return static_cast<Event*>(inventory.get(n));
    }
    
    /// Monte-Carlo simulation step for every Object
    void        step();

};


#endif