// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef SPACE_DISC_H
#define SPACE_DISC_H

#include "space.h"

/// A disc centered at the origin, with variable radius.
/**
 Space `disc` is a disc centered around the origin.
 Forces registered with 'setInteractions' are added, and used to update the
 radius of the Space. How fast the radius changes is set by the value 'mobility'
 in SpaceProp.
 
 Parameters:
     - radius = radius of the disc

 @ingroup SpaceGroup
 
 FJN, Strasbourg 29.01.2017
 */

class SpaceDisc : public Space
{
private:
    
    /// the radius of the disc
    real   radius_;
    
    /// radial force
    mutable real force_;
    
public:
    
    /// constructor
    SpaceDisc(SpaceProp const*);

    /// change dimensions
    void        resize(Glossary& opt);
 
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;

    /// a random position inside the volume
    Vector      randomPlace() const { return Vector::randB(radius_); }
    
    /// direct normal direction calculation
    Vector      normalToEdge(Vector const& pos) const { return normalize(pos); }

    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;

    
    /// add interactions to a Meca
    void        setInteractions(Meca &, FiberSet const&) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
    
    ///    the step function can change the radius
    void        step();

    /// OpenGL display function; returns true if successful
    bool        draw() const;
    
    /// write to file
    void        write(Outputter&) const;

    /// get dimensions from array `len`
    void        setLengths(const real len[8]);
    
    /// read from file
    void        read(Inputter&, Simul&, ObjectTag);

};

#endif

