
// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_squareOffset.h"
#include "exceptions.h"
#include "mecapoint.h"
#include "iowrapper.h"
#include "glossary.h"
#include "meca.h"

SpaceSquareOffset::SpaceSquareOffset(SpaceProp const* p):
SpaceSquare(p)
{
    offset = Vector(0,0,0);
}

void SpaceSquareOffset::resize(Glossary& opt)
{
    SpaceSquare::resize(opt);
    real off[DIM];
    for ( int d = 0; d < DIM; ++d )
    {
        opt.set(off[d], "offset", d);
        if ( off[d] < 0 )
            throw InvalidParameter("squareOffset:offset must be >= 0");
    }
    offset = Vector(off);
}

void SpaceSquareOffset::boundaries(Vector& inf, Vector& sup) const
{
    SpaceSquare::boundaries(inf, sup);
    inf+=offset;
    sup+=offset;
}

bool SpaceSquareOffset::inside(Vector const& w) const
{
    return SpaceSquare::inside(w-offset);
}

Vector SpaceSquareOffset::randomPlace() const
{
    return Space::randomPlace();
}

Vector SpaceSquareOffset::project(Vector const& w) const
{
    return SpaceSquare::project(w-offset)+offset;
}

void SpaceSquareOffset::setInteraction(const real pos[], Mecapoint const& pe, Meca & meca, real stiff, const real dim[], const real off[])
{
    bool in = true;
    
    index_t inx = DIM * pe.matIndex();
    
    real pos_shifted[DIM];
    
    for ( int d = 0; d < DIM; ++d )
    {
        pos_shifted[d] = pos[d]-off[d];
        assert_true( dim[d] >= 0 );
        if ( fabs(pos_shifted[d]) > dim[d] )
        {
            meca.mC(inx+d, inx+d) -= stiff;
            if (pos_shifted[d]>=0)
                meca.base(inx+d)      += stiff * (dim[d]+off[d]);
            else
                meca.base(inx+d)      += stiff * (-dim[d]+off[d]);
            in = false;
        }
    }

    if ( in )
    {
        // find the dimensionality 'dip' corresponding to the closest face
        int  dip = 0;
        
        real l = dim[0] - fabs(pos_shifted[0]);
#if ( DIM > 1 )
        real u = dim[1] - fabs(pos_shifted[1]);
        if ( u < l ) { dip = 1; l = u; };
#endif
#if ( DIM > 2 )
        u = dim[2] - fabs(pos_shifted[2]);
        if ( u < l )  dip = 2;
#endif
        meca.mC(inx+dip, inx+dip) -= stiff;
        if (pos_shifted[dip]>=0)
            meca.base(inx+dip)      += stiff * (dim[dip]+off[dip]);
        else
            meca.base(inx+dip)      += stiff * (-dim[dip]+off[dip]);
    }
}

void SpaceSquareOffset::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    setInteraction(pos, pe, meca, stiff, length_,offset);
}


//------------------------------------------------------------------------------

void SpaceSquareOffset::write(Outputter& out) const
{
    out.put_characters("squareOffset", 16);
    out.writeUInt16(7);
    out.writeFloat(length_[0]);
    out.writeFloat(length_[1]);
    out.writeFloat(length_[2]);
    out.writeFloat(offset[0]);
    out.writeFloat(offset[1]);
    out.writeFloat(offset[2]);
    out.writeFloat(0.f);
}

void SpaceSquareOffset::setLengths(const real len[])
{
    SpaceSquare::setLengths(len);
    
    real off[3];
    off[0] = len[3];
    off[1] = len[4];
    off[2] = len[5];
    offset = off;
    
}

void SpaceSquareOffset::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "squareOffset");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

using namespace gle;

bool SpaceSquareOffset::draw() const
{
    const real X = length_[0] + offset[0];
    const real Y = length_[1] + offset[1];
    const real Z = ( DIM > 2 ) ? length_[2] + offset[2] : 0;
    
    const real Xm = -length_[0] + offset[0];
    const real Ym = -length_[1] + offset[1];
    const real Zm = ( DIM > 2 ) ? -length_[2]+ offset[2] : 0;
  
    glPushAttrib(GL_ENABLE_BIT);

#if ( DIM > 2 )

    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, Zm );
    gleVertex(  X,  Y,  Z );
    gleVertex(  X, Ym, Zm );
    gleVertex(  X, Ym,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex( Xm, Ym, Zm );
    gleVertex( Xm, Ym,  Z );
    gleVertex( Xm,  Y, Zm );
    gleVertex( Xm,  Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, Zm );
    gleVertex( Xm,  Y, Zm );
    gleVertex(  X,  Y,  Z );
    gleVertex( Xm,  Y,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X, Ym,  Z );
    gleVertex( Xm, Ym,  Z );
    gleVertex(  X, Ym, Zm );
    gleVertex( Xm, Ym, Zm );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y,  Z );
    gleVertex( Xm,  Y,  Z );
    gleVertex(  X, Ym,  Z );
    gleVertex( Xm, Ym,  Z );
    glEnd();
    
    glBegin(GL_TRIANGLE_STRIP);
    gleVertex(  X,  Y, Zm );
    gleVertex(  X, Ym, Zm );
    gleVertex( Xm,  Y, Zm );
    gleVertex( Xm, Ym, Zm );
    glEnd();

    glDisable(GL_LIGHTING);
    glLineWidth(0.5);
    glColor3f(1,1,1);
    
    glBegin(GL_LINE_LOOP);
    gleVertex(  X,  Y, Z );
    gleVertex(  X, Ym, Z );
    gleVertex( Xm, Ym, Z );
    gleVertex( Xm,  Y, Z );
    glEnd();
    
    glBegin(GL_LINES);
    gleVertex(  X,  Y, Zm );
    gleVertex(  X,  Y,  Z );
    gleVertex(  X, Ym, Zm );
    gleVertex(  X, Ym,  Z );
    gleVertex( Xm, Ym, Zm );
    gleVertex( Xm, Ym,  Z );
    gleVertex( Xm,  Y, Zm );
    gleVertex( Xm,  Y,  Z );
    glEnd();

#endif
  
    glDisable(GL_LIGHTING);
    glBegin(GL_LINE_LOOP);
    gleVertex(  X,  Y, Zm );
    gleVertex(  X, Ym, Zm );
    gleVertex( Xm, Ym, Zm );
    gleVertex( Xm,  Y, Zm );
    glEnd();

    glPopAttrib();
    return true;
}

#else

bool SpaceSquareOffset::draw() const
{
    return false;
}

#endif
