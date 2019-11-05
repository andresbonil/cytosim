// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "space_strip.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"


SpaceStrip::SpaceStrip(SpaceProp const* p)
: Space(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("strip is not usable in 1D");
    for ( int d = 0; d < 3; ++d )
        length_[d] = 0;
}


void SpaceStrip::resize(Glossary& opt)
{
    for ( int d = 0; d < DIM; ++d )
    {
        real len = length_[d];
        if ( opt.set(len, "length", d) )
            len *= 0.5;
        if ( len < 0 )
            throw InvalidParameter("square:length[] must be >= 0");
        length_[d] = len;
    }
    if ( length_[DIM-1] <= 0 )
    {
        real rad = 0;
        if ( opt.set(rad, "radius") )
            length_[DIM-1] = rad;
        if ( length_[DIM-1] <= 0 )
            throw InvalidParameter("strip:length[DIM-1] must be > 0");
    }
}


Modulo * SpaceStrip::makeModulo() const
{
    Modulo * mod = new Modulo();
    for ( int d = 0; d < DIM-1; ++d )
        mod->enable(d, length_[d]);
    return mod;
}


void SpaceStrip::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_[0],-length_[1],-length_[2]);
    sup.set( length_[0], length_[1], length_[2]);
}

//------------------------------------------------------------------------------
#pragma mark -


#if ( DIM == 1 )

real SpaceStrip::volume() const
{
    return 2.0 * length_[0];
}

bool  SpaceStrip::inside(Vector const& point) const
{
    if ( point[0] >  length_[0] ) return false;
    if ( point[0] < -length_[0] ) return false;
    return true;
}


Vector SpaceStrip::project(Vector const& pos) const
{
    return Vector(std::copysign(length_[0], pos.XX));
}

#endif


//------------------------------------------------------------------------------

#if ( DIM == 2 )

real SpaceStrip::volume() const
{
    return 4.0 * length_[0] * length_[1];
}

bool  SpaceStrip::inside(Vector const& point) const
{
    if ( point[1] >  length_[1] ) return false;
    if ( point[1] < -length_[1] ) return false;
    return true;
}


Vector SpaceStrip::project(Vector const& pos) const
{
    return Vector(pos.XX, std::copysign(length_[1], pos.YY));
}

#endif

//------------------------------------------------------------------------------

#if ( DIM >= 3 )

real SpaceStrip::volume() const
{
    return 8.0 * length_[0] * length_[1] * length_[2];
}

bool  SpaceStrip::inside(Vector const& point) const
{
    if ( point[2] >  length_[2] ) return false;
    if ( point[2] < -length_[2] ) return false;
    return true;
}

Vector SpaceStrip::project(Vector const& pos) const
{
    return Vector(pos.XX, pos.YY, std::copysign(length_[2], pos.ZZ));
}

#endif

//------------------------------------------------------------------------------

void SpaceStrip::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    index_t inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;

#if ( DIM == 2 )
    meca.base(inx) += stiff * std::copysign(length_[1], pos.YY);
#elif ( DIM > 2 )
    meca.base(inx) += stiff * std::copysign(length_[2], pos.ZZ);
#endif
}


void SpaceStrip::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    index_t inx = DIM-1 + DIM * pe.matIndex();
    
    meca.mC(inx, inx) -= stiff;

#if ( DIM == 2 )
    meca.base(inx) += stiff * std::copysign(length_[1]-rad, pos.YY);
#elif ( DIM > 2 )
    meca.base(inx) += stiff * std::copysign(length_[2]-rad, pos.ZZ);
#endif
}

//------------------------------------------------------------------------------

void SpaceStrip::write(Outputter& out) const
{
    out.put_characters("strip", 16);
    out.writeUInt16(4);
    out.writeFloat(length_[0]);
    out.writeFloat(length_[1]);
    out.writeFloat(length_[2]);
    out.writeFloat(0.f);
}


void SpaceStrip::setLengths(const real len[])
{
    length_[0] = len[0];
    length_[1] = len[1];
    length_[2] = len[2];
}

void SpaceStrip::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "strip");
    setLengths(len);
}


//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#pragma mark -

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"
using namespace gle;

bool SpaceStrip::draw() const
{
    const real X = length_[0];
    const real Y = ( DIM > 1 ) ? length_[1] : 1;
    const real Z = ( DIM > 2 ) ? length_[2] : 0;
    
#if ( DIM > 2 )
    // draw faces:
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0, 0, 1);
    gleVertex( -X,  Y, -Z );
    gleVertex(  X,  Y, -Z );
    gleVertex( -X, -Y, -Z );
    gleVertex(  X, -Y, -Z );
    glEnd();
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3f(0, 0, -1);
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    glEnd();
    // draw outline:
    glBegin(GL_LINE_STRIP);
    gleVertex( -X,  Y, -Z );
    gleVertex(  X,  Y, -Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X, -Y, -Z );
    gleVertex( -X,  Y, -Z );
    glEnd();
    glBegin(GL_LINE_STRIP);
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex( -X,  Y, Z );
    glEnd();
#else
    glBegin(GL_LINES);
    gleVertex( -X,  Y, Z );
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex( -X, -Y, Z );
    glEnd();
#endif

    // draw outline:
    glLineStipple(1, 0x000F);
    glEnable(GL_LINE_STIPPLE);
    glBegin(GL_LINES);
#if ( DIM > 2 )
    gleVertex(  X,  Y,  Z );
    gleVertex(  X,  Y, -Z );
    gleVertex(  X, -Y,  Z );
    gleVertex(  X, -Y, -Z );
    gleVertex( -X,  Y,  Z );
    gleVertex( -X,  Y, -Z );
    gleVertex( -X, -Y,  Z );
    gleVertex( -X, -Y, -Z );
#else
    gleVertex(  X,  Y, Z );
    gleVertex(  X, -Y, Z );
    gleVertex( -X,  Y, Z );
    gleVertex( -X, -Y, Z );
#endif
    glEnd();
    glDisable(GL_LINE_STIPPLE);

    return true;
}

#else

bool SpaceStrip::draw() const
{
    return false;
}

#endif

