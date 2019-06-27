// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "vector3.h"
#include "quaternion.h"
#include "gle_color.h"
#include "iowrapper.h"

/// something
class Thing
{
    
private:
    
    /// display color
    gle_color     color;
    
    /// position
    Vector3       pos;
    
public:
    
    ///constructor
    Thing();
    
    ///desctructor
    ~Thing();
    
    ///reset particle
    void reset();
    
    ///perform one simulation step
    void step();
    
    /// draw using OpenGL
    void        draw() const;
    
    ///record state to FILE, the current time is provided
    void write(FILE*, real time) const;

};
