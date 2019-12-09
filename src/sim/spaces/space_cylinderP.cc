// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "space_cylinderP.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "mecapoint.h"
#include "glossary.h"
#include "meca.h"


SpaceCylinderP::SpaceCylinderP(SpaceProp const* p)
: Space(p)
{
    if ( DIM < 3 )
        throw InvalidParameter("cylinderP is only valid in 3D: use strip instead");
    length_ = 0;
    radius_ = 0;
}

void SpaceCylinderP::resize(Glossary& opt)
{
    real len = length_, rad = radius_;
    
    if ( opt.set(rad, "diameter") )
        rad *= 0.5;
    else opt.set(rad, "radius");

    if ( opt.set(len, "length") )
        len *= 0.5;

    if ( rad < 0 )
        throw InvalidParameter("cylinderP:radius must be >= 0");

    if ( len <= 0 )
        throw InvalidParameter("cylinderP:length must be > 0");
    
    length_ = len;
    radius_ = rad;
}


Modulo * SpaceCylinderP::makeModulo() const
{
    Modulo * mod = new Modulo();
    mod->enable(0, length_);
    return mod;
}


void SpaceCylinderP::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-length_,-radius_,-radius_);
    sup.set( length_, radius_, radius_);
}


real SpaceCylinderP::volume() const
{
    return 2 * M_PI * length_ * radius_ * radius_;
}


bool SpaceCylinderP::inside(Vector const& w) const
{
#if ( DIM > 2 )
    const real RT = w.YY * w.YY + w.ZZ * w.ZZ;
    return ( RT <= radius_ * radius_ );
#elif ( DIM > 1 )
    return ( fabs(w.YY) <= radius_ );
#else
    return false;
#endif
}


bool SpaceCylinderP::allInside(Vector const& w, const real rad ) const
{
    assert_true( rad >= 0 );
#if ( DIM > 2 )
    const real RT = w.YY * w.YY + w.ZZ * w.ZZ;
    return ( RT <= square(radius_-rad) );
#elif ( DIM > 1 )
    return ( fabs(w.YY) <= radius_-rad );
#else
    return false;
#endif
}


Vector SpaceCylinderP::randomPlace() const
{
#if ( DIM >= 3 )
    const Vector2 V = Vector2::randB(radius_);
    return Vector(length_*RNG.sreal(), V.XX, V.YY);
#elif ( DIM > 1 )
    return Vector(length_*RNG.sreal(), radius_*RNG.sreal());
#else
    return Vector(length_*RNG.sreal());
#endif
}

//------------------------------------------------------------------------------
Vector SpaceCylinderP::project(Vector const& w) const
{
    Vector p;
    p.XX = w.XX;
    
#if ( DIM > 2 )
    real n = w.normYZ();
    if ( n > REAL_EPSILON )
    {
        p.YY = w.YY * ( radius_ / n );
        p.ZZ = w.ZZ * ( radius_ / n );
    }
    else
    {
        const Vector2 V = Vector2::randU();
        p.YY = radius_ * V.XX;
        p.ZZ = radius_ * V.YY;
    }
#endif
    return p;
}

//------------------------------------------------------------------------------

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    meca.addCylinderClampX(pe, radius_, stiff);
}

/**
 This applies forces towards the cylindrical surface only
 */
void SpaceCylinderP::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    real eRadius = radius_ - rad;
    if ( eRadius < 0 ) eRadius = 0;
    
    meca.addCylinderClampX(pe, eRadius, stiff);
}

//------------------------------------------------------------------------------

void SpaceCylinderP::write(Outputter& out) const
{
    out.put_characters("cylinderP", 16);
    out.writeUInt16(2);
    out.writeFloat(length_);
    out.writeFloat(radius_);
}


void SpaceCylinderP::setLengths(const real len[])
{
    length_ = len[0];
    radius_ = len[1];
}


void SpaceCylinderP::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "cylinderP");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY
#include "opengl.h"
#include "gle.h"

bool SpaceCylinderP::draw() const
{
#if ( DIM > 2 )

    const size_t fin = 512;
    GLfloat c[fin+1], s[fin+1];
    gle::circle(fin, c, s, 1);

    GLfloat L = (GLfloat)length_;
    GLfloat R = (GLfloat)radius_;

    glBegin(GL_TRIANGLE_STRIP);
    for ( size_t n = 0; n <= fin; ++n )
    {
        glNormal3f( 0, c[n], s[n] );
        glVertex3f( +L, R*c[n], R*s[n] );
        glVertex3f( -L, R*c[n], R*s[n] );
    }
    glEnd();
    
    if ( 1 )
    {
        //draw dotted-rings to indicate periodicity
        glLineStipple(1, 0x000F);
        glEnable(GL_LINE_STIPPLE);
        glPushMatrix();
        glTranslatef(L, 0, 0);
        glScalef(R, R, R);
        glRotated(90, 0, 1, 0);
        gle::gleCircle();
        glTranslatef(0, 0, -2*L/R);
        glRotated(180, 0, 1, 0);
        gle::gleCircle();
        glPopMatrix();
        glDisable(GL_LINE_STIPPLE);
    }

#endif
    return true;
}

#else

bool SpaceCylinderP::draw() const
{
    return false;
}

#endif

