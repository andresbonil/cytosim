#include "dim.h"
#include "space_sphereCrown.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "random.h"
#include "meca.h"

SpaceSphereCrown::SpaceSphereCrown(SpaceProp const* p)
: Space(p), inner_radius_(0),outer_radius_(0),inner_radiusSqr_(0),outer_radiusSqr_(0)
{
}

void SpaceSphereCrown::resize(Glossary& opt)
{
    
    opt.set(inner_radius_, "inner_radius");
    opt.set(outer_radius_, "outer_radius");
    
    if ( inner_radius_ < 0||outer_radius_ <0 )
        throw InvalidParameter(prop->name()+":outer_radius and inner_radius must be >= 0");
    if ( inner_radius_ > outer_radius_ )
        throw InvalidParameter(prop->name()+":outer_radius must be bigger than inner_radius");
    update();
}

void SpaceSphereCrown::boundaries(Vector& inf, Vector& sup) const
{
    inf.set(-outer_radius_,-outer_radius_,-outer_radius_);
    sup.set( outer_radius_, outer_radius_, outer_radius_);
}

real SpaceSphereCrown::volume() const
{
#if ( DIM == 1 )
    return 2 * (outer_radius_ - inner_radius_);
#elif ( DIM == 2 )
    return M_PI * (square(outer_radius_)-square(inner_radius_));
#else
    return 4*M_PI/3.0 * (cube(outer_radius_)-cube(inner_radius_));
#endif
}

bool SpaceSphereCrown::inside(Vector const& pos) const
{
    real dist2centerSqr = pos.normSqr();
    return  (dist2centerSqr>=inner_radiusSqr_) && (dist2centerSqr<=outer_radiusSqr_);
}

Vector SpaceSphereCrown::randomPlace() const
{
    Vector pos;
    do {
        pos = Vector::randB(outer_radius_);
    }
    while (!inside(pos));
    return pos;
}

Vector SpaceSphereCrown::project(Vector const& pos) const
{
    real dist = pos.norm();
    
    // We stablish the average radius as the threshold
    real thresh = (outer_radius_ + inner_radius_)/2;
    
    if (dist<=thresh)
    {
        return pos * ( inner_radius_ / dist );
    }
    else
    {
        return pos * ( outer_radius_ / dist );
    }
}

void SpaceSphereCrown::setInteraction(Vector const& pos, Mecapoint const& pe, Meca & meca, real stiff) const
{
    real dist = pos.norm();

    // We stablish the average radius as the threshold
    real thresh = (outer_radius_ + inner_radius_)/2;
    
    if (dist<=thresh)
        meca.addSphereClamp(pos, pe, Vector(0,0,0), inner_radius_, stiff);
    else
    {
        meca.addSphereClamp(pos, pe, Vector(0,0,0), outer_radius_, stiff);
    }
}

void SpaceSphereCrown::setInteraction(Vector const& pos, Mecapoint const& pe, real rad, Meca & meca, real stiff) const
{
    real dist = pos.norm();

    // We stablish the average radius as the threshold
    real thresh = (outer_radius_ + inner_radius_)/2;
    
    if (dist<=thresh)
        meca.addSphereClamp(pos, pe, Vector(0,0,0), rad-inner_radius_, stiff);
    else
    {
        meca.addSphereClamp(pos, pe, Vector(0,0,0), outer_radius_-rad, stiff);
    }
}

//------------------------------------------------------------------------------

void SpaceSphereCrown::write(Outputter& out) const
{
    out.put_characters("sphereCrown", 16);
    out.writeUInt16(2);
    out.writeFloat(inner_radius_);
    out.writeFloat(outer_radius_);
    out.writeFloat(0.f);
}


void SpaceSphereCrown::setLengths(const real len[])
{
    inner_radius_ = len[0];
    outer_radius_ = len[1];
    update();
}


void SpaceSphereCrown::read(Inputter& in, Simul&, ObjectTag)
{
    real len[8] = { 0 };
    read_data(in, len, "sphereCrown");
    setLengths(len);
}

//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------
#ifdef DISPLAY
#include "gle.h"

bool SpaceSphereCrown::draw() const
{
    
    real radi [2] = {inner_radius_,outer_radius_};
    real radius_;
    for (int i=0; i<2; i++) {
        radius_ = radi[i];
    
#if ( DIM <= 2 )
    
    //number of sections in the quarter-circle
    constexpr size_t fin = ((DIM==2) ? 32 : 8) * gle::finesse;
    
    GLfloat cir[2*fin+2];
    gle::circle(fin, cir, (GLfloat)radius_);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(2, GL_FLOAT, 0, cir);
    glDrawArrays(GL_LINE_STRIP, 0, fin+1);
    glDisableClientState(GL_VERTEX_ARRAY);

#else
    
    GLfloat R = (GLfloat)radius_;
    glPushMatrix();
    glScalef(R, R, R);
    gle::gleSphere8B();
    gle::gleThreeBands(128);
    glPopMatrix();
    
#endif
    }
    return true;
}

#else

bool SpaceSphereCrown::draw() const
{
    return false;
}

#endif
