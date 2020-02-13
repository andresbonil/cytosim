// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_DYNAMIC_ELLIPSE_PROP_H
#define SPACE_DYNAMIC_ELLIPSE_PROP_H

#include "space_prop.h"
#include "property.h"

class SpaceDynamicEllipse;
class SpaceProp;
class Space;


class SpaceDynamicEllipseProp : public SpaceProp 
{
    
    friend class SpaceDynamicEllipse;
    
public:
    
    // tension of the ellipse
    real    tension;
    
    // volume of the ellipse (mutable because changed by const method)
    mutable real    volume;
	
public:
    
	
    /// constructor
    SpaceDynamicEllipseProp(const std::string& n) : SpaceProp(n)  { clear(); }
    
    /// destructor
    ~SpaceDynamicEllipseProp() { }
    
    /// create a new, uninitialized, Space
    Space * newSpace() const;

	/// create a new Space according to specifications
	Space * newSpace(Glossary&) const;
        
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
	
    /// check and derive more parameters
    void complete(Simul const&);
	
    /// write all values
    void write_values(std::ostream&) const;
    
	
};

#endif

