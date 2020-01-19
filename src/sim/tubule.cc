// Cytosim was created by Francois Nedelec.
// Copyright Cambridge University, 2019

#include "dim.h"
#include "tubule.h"
#include "fiber.h"
#include "simul.h"
#include "meca.h"
#include "exceptions.h"
#include "messages.h"


Tubule::Tubule(TubuleProp * p) : prop(p)
{
    reset();
}


Tubule::~Tubule()
{
    prop = nullptr;
}


void Tubule::reset()
{
    bone_ = nullptr;
    for ( size_t i = 0; i < FILM; ++i )
        fil_[i] = nullptr;
}


void Tubule::setFamily(Fiber const* fam)
{
    // wrap array values for convenience
    for ( size_t i = NFIL; i < FILM; ++i )
        fil_[i] = fil_[i-NFIL];

#if FIBER_HAS_FAMILY
    for ( size_t i = 0; i < NFIL; ++i )
    {
        fil_[i]->family_ = fam;
        fil_[i]->sister_ = fil_[(i+NFIL-1)%NFIL];
        fil_[i]->brother_ = fil_[(i+1)%NFIL];
    }
    if ( bone_ )
        bone_->family_ = fam;
#else
    LOG_ONCE("WARNING: Tubule requires FIBER_HAS_FAMILY\n");
#endif
}


ObjectList Tubule::build(Glossary& opt, Simul& sim)
{
    ObjectList res;
    
    Property * fip = sim.properties.find("fiber", prop->bone_type);
    if ( fip )
    {
        bone_ = static_cast<FiberProp*>(fip)->newFiber(opt);
        res.push_back(bone_);
    }
    
    for ( size_t i = 0; i < NFIL; ++i )
    {
        ObjectList objs = sim.fibers.newObjects(prop->fiber_type, opt);
        res.append(objs);
        fil_[i] = Fiber::toFiber(objs[0]);
    }
    
#if ( DIM >= 3 )
    Vector E(0,tube_radius,0), F(0,0,tube_radius);
    {
        // find average direction
        Vector dir(0,0,0);
        for ( size_t i = 0; i < NFIL; ++i )
            dir += fil_[i]->diffPoints(0);
        dir.normalize();
        dir.orthonormal(E, F);
        E *= tube_radius;
        F *= tube_radius;
    }
    // adjust protofilaments to form a tube:
    real a = M_PI * RNG.sreal();
    real da = 2 * M_PI / NFIL;
    for ( size_t i = 0; i < NFIL; ++i )
    {
        fil_[i]->translate(cos(a)*E+sin(a)*F);
        a += da;
    }
#endif
    
    // set as left-handed helix:
    assert_true(signature());
    for ( size_t i = 0; i < NFIL; ++i )
    {
        Buddy::connect(fil_[i]);
        fil_[i]->setOrigin(i*(-0.012/NFIL));
    }

    setFamily(bone_?bone_:fil_[0]);
    
    return res;
}


void Tubule::goodbye(Buddy * b)
{
    std::cerr << "ERROR: Tubule's filaments cannot be deleted\n";
    if ( b )
    {
        for ( size_t i = 0; i < NFIL+2; ++i )
            if ( fil_[i] == b )
            {
                fil_[i] = nullptr;
                return;
            }
        std::cerr << "Error: unknown cytosim buddy" << b << '\n';
    }
}


void Tubule::step(Simul&)
{
}

    
/**
This uses only addSideLink() with appropriate directions
Cambridge, 12.2019
*/
void Tubule::setInteractionsB(Meca& meca)
{
#if ( DIM >= 3 )
    const real stiff = prop->stiffness[0];
    const real ang = M_PI / NFIL;
    const real len = 2 * tube_radius * sin(ang);  // distance between protofilaments
    const real c = cos(ang), s = sin(ang);
    
    assert_true(fil_[0]);
    const size_t end = fil_[0]->nbPoints() - 1;
    
    Rotation mat(0,1);

    for ( size_t i = 0; i < end; ++i )
    {
        // get centerline
        Vector cen(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            cen += fil_[n]->posPoint(i);
        cen /= NFIL;
        
        // get average direction of the Tubule at this location:
        Vector dir(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            dir += fil_[n]->diffPoints(i);
        dir.normalize();
        mat = Rotation::rotationAroundAxis(dir, c, s);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            Vector arm = mat.vecmul(( cen - fil_[n]->posPoint(i) ).normalized(len));
            meca.addSideLink3D(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), arm, stiff);
        }
    }

    // get centerline
    Vector cen(0,0,0);
    for ( size_t n = 0; n < NFIL; ++n )
        cen += fil_[n]->posPoint(end);
    cen /= NFIL;
    
    for ( size_t n = 0; n < NFIL; ++n )
    {
        Vector arm = mat.vecmul(( cen - fil_[n]->posPoint(end) ).normalized(len));
        meca.addSideLink3D(Interpolation(fil_[n],end-1,end,1), Mecapoint(fil_[n+1],end), arm, stiff);
    }
#endif
}


/**
 This uses a centerline `bone`
 Cambridge, 18.01.2020
 */
void Tubule::setInteractions(Meca& meca)
{
#if ( DIM >= 3 )
    const real stiff = prop->stiffness[0];
    const real ang = M_PI / NFIL;
    const real len = 2 * tube_radius * sin(ang);  // distance between protofilaments
    const real c = cos(ang), s = sin(ang);
    
    assert_true(fil_[0]);
    const size_t end = fil_[0]->nbPoints() - 1;
    
    Rotation mat(0,1);
    real alpha = tube_radius / len;
    Vector cen, dir;
    
    size_t i;
    for ( i = 0; i < end; ++i )
    {
        // get centerline
        cen = bone_->posPoint(i);
        dir = bone_->dirSegment(i);
        mat = Rotation::rotationAroundAxis(dir, c, s);
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            Vector arm = ( cen - fil_[n]->posPoint(i) ).normalized(len);
            meca.addSideLink3D(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), mat.vecmul(arm), stiff);
            meca.addSideLink3D(Interpolation(fil_[n],i,i+1,0), Mecapoint(bone_,i), alpha*cross(dir,arm), stiff);
        }
    }

    // get centerline
    cen = bone_->posPoint(i);
    for ( size_t n = 0; n < NFIL; ++n )
    {
        Vector arm = ( cen - fil_[n]->posPoint(i) ).normalized(len);
        meca.addSideLink3D(Interpolation(fil_[n],i-1,i,1), Mecapoint(fil_[n+1],i), mat.vecmul(arm), stiff);
        meca.addSideLink3D(Interpolation(fil_[n],i-1,i,1), Mecapoint(bone_,i), alpha*cross(dir,arm), stiff);
    }
#endif
}


/**
This uses addSideLink() and addTorque() with appropriate directions
Cambridge, 12.2019
*/
void Tubule::setInteractionsC(Meca& meca)
{
    const real stiffL = prop->stiffness[0];
    const real stiffA = prop->stiffness[1];
    const real ang = M_PI / NFIL;
    const real len = 2 * tube_radius * sin(ang);  // distance between protofilaments
#if ( DIM >= 3 )
    real co = cos(2*ang), si = sin(2*ang);
    
    assert_true(fil_[0]);
    const size_t end = fil_[0]->nbPoints() - 1;
    
    MatrixBlock mat;
    for ( size_t i = 0; i <= end; ++i )
    {
        // get centerline
        Vector cen(0,0,0);
        for ( size_t n = 0; n < NFIL; ++n )
            cen += fil_[n]->posPoint(i);
        cen /= NFIL;
        
        if ( i < end )
        {
            // get average direction of the Tubule at this location:
            Vector dir(0,0,0);
            for ( size_t n = 0; n < NFIL; ++n )
                dir += fil_[n]->diffPoints(i);
            dir.normalize();
            
            // create rotation matrix for torque:
            mat = Meca::torqueMatrix(stiffA, dir, co, si);
            
            for ( size_t n = 0; n < NFIL; ++n )
            {
                Vector arm = (2*cen - fil_[n]->posPoint(i)- fil_[n+1]->posPoint(i)).normalized(len);
                meca.addSideLink3D(Interpolation(fil_[n],i,i+1,0), Mecapoint(fil_[n+1],i), arm, stiffL);
            }
        }
        else
        {
            for ( size_t n = 0; n < NFIL; ++n )
            {
                Vector arm = (2*cen - fil_[n]->posPoint(i) - fil_[n+1]->posPoint(i)).normalized(len);
                meca.addSideLink3D(Interpolation(fil_[n],i-1,i,1), Mecapoint(fil_[n+1],i), arm, stiffL);
            }
        }
        
        for ( size_t n = 0; n < NFIL; ++n )
        {
            meca.addTorque(Mecapoint(fil_[n],i), Mecapoint(fil_[n+1],i),
                           Mecapoint(fil_[n+2],i), mat, stiffA);
        }
    }
#endif
}


void Tubule::write(Outputter& out) const
{
    out.writeUInt16(NFIL);
    out.writeSoftNewline();
    for ( size_t i = 0; i < NFIL; ++i )
    {
        out.writeSoftSpace();
        fil_[i]->writeReference(out);
    }
}


void Tubule::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    size_t n = in.readUInt16();
    if ( n != NFIL )
        throw InvalidIO("unexpected number of filaments in Tubule");
    
    ObjectTag g;
    for ( size_t i = 0; i < n; ++i )
    {
        Object * w = sim.readReference(in, g);
        if ( i < NFIL )
            fil_[i] = Fiber::toFiber(w);
    }
    setFamily(fil_[0]);
}

