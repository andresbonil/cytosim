// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SPHERECROWN_H
#define SPACE_SPHERECROWN_H

#include "space.h"

/// sphere crown centered at the origin.
/**
 Space `sphereCrown` is a sphere crown centered around the origin
 
 Parameters:
    - radius_in = inner radius
    - radius_out = outer_radius
    .
 
 @ingroup SpaceGroup
 */

class SpaceSphereCrown : public Space
{
protected:
    
    /// the inner radius
    real  inner_radius_;
    
    /// the outer radius
    real  outer_radius_;
    
    /// square of the inner radius
    real  inner_radiusSqr_;
    
    /// square of the outer radius
    real  outer_radiusSqr_;
    
    /// calculate radiusSqr
    void  update() { inner_radiusSqr_ = square(inner_radius_);outer_radiusSqr_ = square(outer_radius_); }
    
public:
    
    /// constructor
    SpaceSphereCrown(SpaceProp const*);

    /// check number and validity of specified lengths
    void        resize(Glossary& opt);

    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// radius
    real        thickness() const { return outer_radius_ - inner_radius_; }

    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const;
    
    /// direct normal direction calculation
    Vector      normalToEdge(Vector const& pos) const { return normalize(pos); }

    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;

    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
    
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

