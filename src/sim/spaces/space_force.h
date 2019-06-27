// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_FORCE_H
#define SPACE_FORCE_H

#include "space.h"

/// generates a force field
/**
 Space `force` generates a force field that applies to every objects.
 It is sufficient to create a space, for the force to be active.
 This space cannot be used for confinement.
 
 Parameter:
     - force = vector force
     - center = origin of force field
     - stiffness = stiffness
 .
 
 @ingroup SpaceGroup
 */
class SpaceForce : public Space
{
    /// stiffness of interaction
    real        stiffness;
    
    /// center
    Vector      center;
    
    /// force applied in every point
    Vector      force;
    
public:
    
    ///creator
    SpaceForce(SpaceProp const*);
    
    /// change dimensions
    void        resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const { return true; }

    /// true if a sphere (center w, radius) fits in the space, edges included
    bool        allInside(Vector const&, real rad) const { return true; }
    
    /// true if a sphere (center w[], radius) is entirely outside
    bool        allOutside(Vector const&, real rad) const { return false; }
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;

    /// apply force to all objects in Meca
    void        setInteractions(Meca&, FiberSet const&) const;
    
    /// OpenGL display function; returns true if successful
    bool        draw() const;

};

#endif


