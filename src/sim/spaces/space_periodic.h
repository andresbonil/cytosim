// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_PERIODIC_H
#define SPACE_PERIODIC_H

#include "space.h"
#include "modulo.h"

/// a rectangular Space with periodic boundary conditions
/**
 Space `periodic` implements periodic boundary condition in all dimensions.
 The volume has no edge and wraps on itself.

 Parameters:
     - length = total extent in X, Y and Z
     .
 
 To display a periodic Space, use simul:display parameter 'tile'.
 @ingroup SpaceGroup
 */
class SpacePeriodic : public Space
{
    
    /// dimensions
    real   length_[3];
 
public:
    
    /// creator
    SpacePeriodic(SpaceProp const*);

    /// change dimensions
    void        resize(Glossary& opt);
    
    /// initialize Modulo Object
    Modulo *    makeModulo() const;

    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// the volume inside
    real        volume() const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// set `proj` as the point on the edge that is closest to `point`
    Vector      project(Vector const& pos) const;
    
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

