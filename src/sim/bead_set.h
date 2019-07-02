// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include "object_set.h"
#include "bead.h"

class Simul;

/// a list of Bead
class BeadSet : public ObjectSet
{
public:
    
    /// creator
    BeadSet(Simul& s) : ObjectSet(s) {}
    
    //--------------------------
    
    /// identifies the class
    std::string title() const { return "bead"; }
    
    /// create a new property of category `cat` for a class `name`
    Property *  newProperty(const std::string& cat, const std::string& name, Glossary&) const;
    
    /// create objects of class `name`, given the options provided in `opt`
    ObjectList  newObjects(const std::string& name, Glossary& opt);
    
    /// create a new object (used for reading trajectory file)
    Object *    newObjectT(ObjectTag, unsigned);
    
    //--------------------------
    
    /// remove from the list
    void        remove(Object *);
    
    /// first Object
    Bead *      first() const
    {
        return static_cast<Bead*>(nodes.front());
    }
    
    /// first Bead in inventory
    Bead *      firstID() const
    {
        return static_cast<Bead*>(inventory.first());
    }
    
    /// next Bead in inventory
    Bead *      nextID(Bead const* obj) const
    {
        return static_cast<Bead*>(inventory.next(obj));
    }

    /// find object from its Number
    Bead *      findID(ObjectID n) const
    {
        return static_cast<Bead*>(inventory.get(n));
    }
    
    /// modulo the position (periodic boundary conditions)
    void        foldPosition(Modulo const*) const;
    
    /// Monte-Carlo simulation step for every Object
    void        step() {}
};


#endif
