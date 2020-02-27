
// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SPACE_SQUAREOFFSET_H
#define SPACE_SQUAREOFFSET_H

#include "space_square.h"

///a subclass of SpaceSquare that implements a square that is not centered at the origin
/**
 Space `squareOffset` is a a square that is not centered at the origin.
 
 Parameters:
     - length = extent in X, Y and Z
     - offset = vector indicating the displacement from the origin
 .

 @ingroup SpaceGroup
 */
class SpaceSquareOffset : public SpaceSquare
{
    /// offset
    Vector   offset;
    
protected:
    
    /// apply a force directed towards the edge of the Space
    static void setInteraction(const real pos[], Mecapoint const&, Meca &, real stiff, const real dim[], const real off[]);

public:
    
    ///creator
    SpaceSquareOffset(SpaceProp const*);
    
    /// change dimensions
    void        resize(Glossary& opt);
    
    /// return bounding box in `inf` and `sup`
    void        boundaries(Vector& inf, Vector& sup) const;
    
    /// true if the point is inside the Space
    bool        inside(Vector const&) const;
    
    /// a random position inside the volume
    Vector      randomPlace() const;
    
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
    
    /// apply a force directed towards the edge of the Space
    void        setInteraction(Vector const& pos, Mecapoint const&, Meca &, real stiff) const;
    
    /// apply a force directed towards the edge of the Space
//    void        setInteraction(Vector const& pos, Mecapoint const&, real rad, Meca &, real stiff) const;
};

#endif


