// Cytosim was created by Francois Nedelec. Copyright 2021 Cambridge University

#include "dim.h"
#include "assert_macro.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "aster.h"
#include "solid.h"
#include "solid_prop.h"
#include "fiber_prop.h"
#include "pointsonsphere.h"
#include "random_vector.h"
#include "mecapoint.h"
#include "interpolation.h"
#include "glossary.h"
#include "simul.h"
#include "meca.h"


void Aster::step()
{
    // nucleation:
    for ( size_t i = 0; i < asLinks.size(); ++i )
    {
        if ( !fiber(i) &&  RNG.test(prop->fiber_prob) )
        {
            ObjectList objs;
            Simul & sim = simul();
            Vector A = posSolid1(i);
            Vector B = posSolid2(i);
            Fiber * F = makeFiber(objs, sim, A, B-A, prop->fiber_type, prop->fiber_spec);
            grasp(F, i);
            sim.add(objs);
        }
    }
}

#if   ( DIM == 1 )
#    define ADDLINK addLink2
#elif ( DIM == 2 )
#    define ADDLINK addLink3
#else
#    define ADDLINK addLink4
#endif

void Aster::setInteractions(Meca& meca) const
{
    Solid const* sol = solid();
    if ( !sol )
        return;

    for ( size_t n = 0 ; n < asLinks.size(); ++n )
    {
        Fiber * fib = fiber(n);

        if ( fib )
        {
            AsterLink const& link = asLinks[n];
            
            const unsigned off = sol->matIndex() + link.prime_;
            const unsigned pts[] = { off, off+1, off+2, off+3 };

#ifdef BACKWARD_COMPATIBILITY
            if ( link.alt_ > 0 )
            {
                meca.addLink(Mecapoint(sol, link.prime_), fib->exactEnd(prop->focus), prop->stiffness[0]);
                if ( fib->length() > link.len_ )
                {
                    meca.addLink(Mecapoint(sol, link.alt_), fib->interpolate(link.len_, prop->focus), prop->stiffness[1]);
                }
                else
                {
                    FiberEnd tip = ( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
                    // link the opposite end to an interpolation of the two solid-points:
                    real c = fib->length() / link.len_;
                    meca.addLink(fib->exactEnd(tip), Interpolation(sol, link.prime_, link.alt_, c), prop->stiffness[1]);
                }
                continue;
            }
#endif
            if ( link.rank_ == 1 )
                meca.addLink(fib->exactEnd(prop->focus), Mecapoint(sol, link.prime_), prop->stiffness[0]);
            else
                meca.ADDLINK(fib->exactEnd(prop->focus), pts, link.coef1_, prop->stiffness[0]);
            
            
            // make second type of link:
            real len = link.len_;
            
            if ( fib->length() >= len )
            {
                if ( len > 0 )
                    meca.ADDLINK(fib->interpolate(len, prop->focus), pts, link.coef2_, prop->stiffness[1]);
                else
                    meca.ADDLINK(fib->exactEnd(prop->focus), pts, link.coef2_, prop->stiffness[1]);
            }
            else
            {
                // link the opposite fiber end to a new interpolation:
                FiberEnd end = ( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
                real c = fib->length() / len;
                real u = 1.0 - c;
                real coef[4];
                for ( int d = 0; d < 4; ++d )
                    coef[d] = u * link.coef1_[d] + c * link.coef2_[d];
                meca.ADDLINK(fib->exactEnd(end), pts, coef, prop->stiffness[1]);
            }
        }
    }
}


Aster::~Aster()
{
    //Cytosim::log("destroying %c%lu\n", TAG, identity());
    asSolid = nullptr;
    prop = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Build

/**
 @defgroup NewAster How to create an Aster
 @ingroup NewObject
 
 By default the aster creates a radial distribution of fiber,
 and only the radius need to be specified:
 
     new aster NAME
     {
       fibers = INTEGER, FIBER_NAME, FIBER_SPEC
       radius = OUTER_RADIUS, INNER_RADIUS
       ...
     }
 
 The configuration of the Aster can also be customized by specifying
 directly the points on which the fibers are attached:
 
     new NAME
     {
       radius = 1
       position = -3 0 to 3 0
       point1 = 0 0 0, 0.2
       fiber1 = 0 -0.1 0, 0.1 -0.1 0
       fiber2 = 0  0.0 0, 0.1  0.0 0
       fiber3 = 0  0.1 0, 0.1  0.1 0
     }
 
 One can use an existing Solid to build an aster:
 
     new solid core
     {
     ...
     }
 
     new NAME
     {
        radius = 1
        solid = core1
     }
 
 
 The aster 'type' can be:
 - `astral` fiberd are anchored at random positions near the center, pointing outward
 - `radial` fibers are anchored always at the same distance from the center, pointing radially
 - `regular` fibers are anchored regularly over the surface and point radially
 - `angular` where all fibers are restricted within an specified solid angle,
 .

 */
ObjectList Aster::build(Glossary& opt, Simul& sim)
{
    assert_true(prop);
    assert_true(asSolid==nullptr);
    assert_true(nbOrganized()==0);

    ObjectList objs;
    opt.set(asRadius, "radius");
    if ( asRadius <= 0 )
        throw InvalidParameter("aster:radius must be specified and > 0");
    
    if ( opt.has_value("nb_fibers", 0) )
    {
        opt.define("fiber", 2, opt.value("fiber", 1));
        opt.define("fiber", 1, opt.value("fiber", 0));
        opt.define("fiber", 2, opt.value("nb_fibers", 0));
        throw InvalidParameter("please specify `fibers = NUMBER, CLASS, SPEC`");
    }

    size_t origin = makeSolid(objs, sim, opt);
    
    unsigned type = 7 * opt.has_key("fiber1");
    opt.set(type, "type", {{"radial", 0}, {"astral", 1}, {"regular", 2},
                           {"angular", 3}, {"disc", 4}, {"custom", 7}});
    switch( type )
    {
        case 0: build0(objs, opt, sim, origin); break;
        case 1: build1(objs, opt, sim, origin); break;
        case 2: build2(objs, opt, sim, origin); break;
        case 3: build3(objs, opt, sim, origin); break;
        case 4: build4(objs, opt, sim, origin); break;
        case 7: build7(objs, opt, sim, origin); break;
        default:
            throw InvalidParameter("unknown aster:type");
    }
    return objs;
}


Fiber * Aster::makeFiber(ObjectList& res, Simul& sim, const Vector pos, Vector dir,
                         std::string const& fiber_type, std::string const& fos)
{
    real n = dir.normSqr();
    
    if ( n > REAL_EPSILON )
        dir /= std::sqrt(n);
    else
        dir = Vector::randU();
    
    if ( prop->focus == PLUS_END )
        dir = -dir;
    
    if ( fos.empty() )
        throw InvalidParameter("Error: unspecified fiber specs (aster:fibers[2])");
    
    Glossary opt(fos);
    ObjectList objs = sim.fibers.newObjects(fiber_type, opt);
    
    if ( objs.empty() )
        throw InvalidParameter("could not create aster:fiber");

    Fiber * F = Fiber::toFiber(objs[0]);

    if ( !F )
        throw InvalidParameter("unexpected object returned by fibers.newObjects()");

    opt.print_warning(std::cerr, 1, "aster:build\n");

    //std::clog << "new aster:fiber " << pos << " and " << dir << "\n";
    ObjectSet::rotateObjects(objs, Rotation::rotationToVector(dir));
    ObjectSet::translateObjects(objs, asRadius*pos - F->posEnd(prop->focus));
    
    res.append(objs);
    return F;
}


/**
 define new Anchor point for a Fiber between at positions A and B, specified
 in a local reference frame associated with the Aster: (0,0,0) is 'ref' and
 (1,0,0) is 'ref+1' and (0,1,0) is 'ref+2'.
 Dimensions will be scaled by 'asRadius' because that is the distance between
 the Solid points.
 */
size_t Aster::placeAnchor(const Vector A, const Vector B, size_t ref)
{
    AsterLink & link = asLinks.new_val();
    //std::clog << "Aster::placeAnchor(" << asLinks.size() << ")\n";
    link.set(A, B, ref);
    link.len_ *= asRadius;
    //link.print(std::clog);
    return asLinks.size() - 1;
}


void Aster::placeFiber(ObjectList& objs, Simul& sim, const Vector A, const Vector B,
                       size_t ref, std::string const& fiber_type, std::string const& fos)
{
    size_t i = placeAnchor(A, B, ref);
    Fiber * F = makeFiber(objs, sim, A, B-A, fiber_type, fos);
    assert_true( i == nbOrganized() );
    nbOrganized(i+1);
    grasp(F, i);
}


size_t Aster::makeSolid(ObjectList& objs, Simul& sim, Glossary& opt)
{
    Solid * sol = nullptr;
    // find the Solid specified:
    std::string spec;
    if ( opt.set(spec, "solid") )
    {
        SolidProp const* p = static_cast<SolidProp*>(sim.findProperty("solid", spec));
        
        if ( p )
        {
            sol = new Solid(p);
            objs = sol->build(opt, sim);
            objs.push_back(sol);
            //std::clog << "Aster::makeSolid() created solid " << sol->reference() << "\n";
        }
        else
        {
            sol = Solid::toSolid(sim.solids.findObject(spec, "solid"));
            if ( sol )
            {
                // prevent Aster from being moved, so that its position match the Solid
                opt.define("placement", "off");
                //std::clog << "Aster created on solid " << sol->reference() << "\n";
            }
            else
                throw InvalidParameter("unknown aster:solid `"+spec+"'");
        }
    }
    else
        throw InvalidParameter("aster:solid must be specified");
    
    // check that there is at least one point:
    if ( sol->dragCoefficient() < REAL_EPSILON )
#ifdef BACKWARD_COMPATIBILITY
        sol->addSphere(Vector(0,0,0), asRadius);
#else
        throw InvalidParameter("Aster's drag coefficient is null: please specify 'point1=center, RADIUS'");
#endif
    
    // add local coordinate system around the last point:
    size_t ref = sol->addTriad(asRadius);
    sol->fixShape();
    asSolid = sol;
    //asSolid->write(std::clog);
    if ( !Buddy::check(asSolid) )
        Buddy::connect(asSolid);
    return ref;
}


/// fiber's anchor points specified directly
void Aster::build7(ObjectList& objs, Glossary& opt, Simul& sim, size_t ref)
{
    std::string tif, fos;
    opt.set(tif, "fibers");
    opt.set(fos, "fibers", 1);

    Vector A, B;
    size_t cnt = 1;
    std::string var = "fiber1";
    while ( opt.set(A, var) && opt.set(B, var, 1) )
    {
        //std::clog << "direct fiber anchor " << pos1 << " and " << pos2 << "\n";
        std::string str = fos;
        opt.set(str, var, 2);
        placeFiber(objs, sim, A, B, ref, tif, str);
        var = "fiber" + std::to_string(++cnt);
    }
}


/// This is a special case for Yeast's Spindle Pole Bodies
void Aster::build4(ObjectList& objs, Glossary& opt, Simul& sim, size_t ref)
{
    real dis = 0;
    real sep = 0.025; // 25 nm by default, corresponding to Microtubules
    size_t nbf = 7;
    std::string tif, fos;
    opt.set(nbf, "fibers");
    opt.set(tif, "fibers", 1);
    opt.set(fos, "fibers", 2);
    opt.set(dis, "radius", 1);
    dis *= 0.5;
    opt.set(sep, "seed_diameter");
    size_t ouf = 0;
    size_t cnt = 0;
    std::vector<Vector2> pts(nbf, Vector2(0,0));
    do {
        cnt = tossPointsDisc(pts, sep/asRadius, 1024);
    } while ( cnt < nbf && ++ouf < 1024 );
    if ( cnt < nbf )
    {
        std::clog << "warning: aster could only fit " << cnt << " seeds ";
        std::clog << "with aster:seed_diameter = " << sep << '\n';
    }
    //std::clog << "toss(" << nbf << ") placed " << cnt << "\n";
    for ( size_t i = 0; i < cnt; ++i )
    {
        // orient anchors by default along the X-axis:
        real y = pts[i].y();
#if ( DIM == 3 )
        real x = pts[i].XX;
        Vector3 A(-dis, x, y);
        Vector3 B( dis, x, y);
#else
        Vector A(-dis, y, 0);
        Vector B( dis, y, 0);
#endif
        placeFiber(objs, sim, A, B, ref, tif, fos);
    }
}


/**
 For type 'angular' all fibers are restricted within an specified solid angle,
 and their orientation is radial
 initial code by GAELLE LETORT, 14.03.2017
 */
void Aster::build3(ObjectList& objs, Glossary& opt, Simul& sim, size_t ref)
{
    real dis = 0;
    size_t nbf = 7;
    std::string tif, fos;
    opt.set(nbf, "fibers");
    opt.set(tif, "fibers", 1);
    opt.set(fos, "fibers", 2);
    opt.set(dis, "radius", 1);
    dis /= asRadius;

    real cap = 1.0, angle = M_PI;
    // either 'angle' or 'cap' can be specified:
    if ( opt.set(angle, "aster_angle") )
        cap = 1.0 - std::cos(angle);
    else
        opt.set(cap, "aster_cap" );

    std::vector<Vector> pts(nbf, Vector(0,0,0));
    size_t cnt = nbf;

#if ( DIM == 1 )
    // No effect of angle in 1D:
    for ( size_t i = 0; i < nbf; ++i )
        pts[i] = Vector1(2*(i&2)-1);
#elif ( DIM == 2 )
    real delta = 2 * angle / real(nbf);
    // points are evenly distributed in [-aster_angle, aster_angle]
    real ang = -angle;
    for ( size_t i = 0; i < nbf; ++i )
    {
        pts[i] = Vector(std::cos(ang), std::sin(ang));
        ang += delta;
    }
#else
    // distribute points randomly over a portion of the unit sphere:
    size_t ouf = 0;
    real sep, sep0 = std::sqrt( 2 * M_PI * cap / nbf );
    do {
        // we decrease gradually the separation, to reach a good solution...
        sep = 512 * sep0 / real(ouf+512);
        cnt = tossPointsCap(pts, cap, sep, 1024);
        //std::clog << "tossCap(" << nbf << ") placed " << cnt << " with sep = " << sep << "\n";
    } while ( cnt < nbf && ++ouf < 1024 );
    if ( cnt < nbf )
        std::clog << "warning: aster could only fit " << cnt << " seeds\n";
    //std::clog << "tossCap(" << nbf << ") placed " << cnt << " with sep = " << sep << "\n";
#endif
    
    for ( size_t i = 0; i < cnt; ++i )
    {
        Vector B = pts[i];
        placeFiber(objs, sim, dis*B, B, ref, tif, fos);
    }
}

/**
 For type 'regular', fibers are regularly distributed on the surface,
 */
void Aster::build2(ObjectList& objs, Glossary& opt, Simul& sim, size_t ref)
{
    real dis = 0;
    size_t nbf = 7;
    std::string tif, fos;
    opt.set(nbf, "fibers");
    opt.set(tif, "fibers", 1);
    opt.set(fos, "fibers", 2);
    opt.set(dis, "radius", 1);
    dis /= asRadius;
#if ( DIM == 1 )
    Vector A(0, 0, 0);
    Vector B(dis, 0, 0);
    for ( size_t i = 0; i < nbf; i += 2 )
    {
        real S = ( i & 1 ) ? -1 : 1;
        placeFiber(objs, sim, A, S*B, ref, tif, fos);
    }
#elif ( DIM == 2 )
    real a = 0;
    real delta = 2 * M_PI / real(nbf);
    for ( size_t i = 0; i < nbf; ++i, a+=delta )
    {
        Vector B(std::cos(a), std::sin(a));
        placeFiber(objs, sim, dis*B, B, ref, tif, fos);
    }
#else
    //we use SphericalCode to distribute points 'equally' on the sphere
    PointsOnSphere code(nbf);
    Vector B(0, 0, 0);
    for ( size_t i = 0; i < nbf; ++i )
    {
        code.copyPoint(B, i);
        placeFiber(objs, sim, dis*B, B, ref, tif, fos);
    }
#endif
}


/**
 For type 'astral' we put fibers randomly on the surface,
 with a constrain based on the scalar product: position*direction > 0
 */
void Aster::build1(ObjectList& objs, Glossary& opt, Simul& sim, size_t ref)
{
    real dis = 0;
    size_t nbf = 7;
    std::string tif, fos;
    opt.set(nbf, "fibers");
    opt.set(tif, "fibers", 1);
    opt.set(fos, "fibers", 2);
    opt.set(dis, "radius", 1);
    dis /= asRadius;
    for ( size_t i = 0; i < nbf; ++i )
    {
        Vector P = Vector::randB();
        Vector D = Vector::randU();
        while ( dot(D, P) < 0 )
            D = Vector::randU();
        placeFiber(objs, sim, P-dis*D, P, ref, tif, fos);
    }
}


/**
 For type 'radial' we put fibers randomly on the surface, and set their
 direction as purely radial. We require a separation of 25 nm by default,
 corresponding to Microtubule's size.
 */
void Aster::build0(ObjectList& objs, Glossary& opt, Simul& sim, size_t ref)
{
    real dis = 0;
    real sep = 0.025; // 25 nm by default, corresponding to Microtubules
    size_t nbf = 7;
    std::string tif, fos;
    opt.set(nbf, "fibers");
    opt.set(tif, "fibers", 1);
    opt.set(fos, "fibers", 2);
    opt.set(dis, "radius", 1);
    opt.set(sep, "seed_diameter");
    dis /= asRadius;
    size_t ouf = 0;
    size_t cnt = 0;
    std::vector<Vector> pts(nbf, Vector(0,0,0));
    do {
        cnt = tossPointsSphere(pts, sep/asRadius, 1024);
    } while ( cnt < nbf && ++ouf < 1024 );
    if ( cnt < nbf )
    {
        std::clog << "warning: aster could only fit " << cnt << " seeds ";
        std::clog << "with aster:seed_diameter = " << sep << '\n';
    }
    //std::clog << "toss(" << nbf << ") placed " << cnt << "\n";
    for ( size_t i = 0; i < cnt; ++i )
    {
        Vector B = pts[i];
        placeFiber(objs, sim, dis*B, B, ref, tif, fos);
    }
}


//------------------------------------------------------------------------------
#pragma mark - I/O

void Aster::write(Outputter& out) const
{
    out.writeUInt16(nbOrganized()+1);
    out.writeSoftNewline();
    Object::writeReference(out, asSolid);
    for ( size_t i = 0; i < nbOrganized(); ++i )
    {
        out.writeSoftSpace();
        Object::writeReference(out, organized(i));
    }
    //Organizer::write(out);
    
    out.writeSoftNewline();
    out.writeUInt16(asLinks.size());
    for ( size_t ii = 0; ii < asLinks.size(); ++ii )
    {
        out.writeSoftNewline();
        asLinks[ii].write(out);
    }
}


void Aster::read(Inputter& in, Simul& sim, ObjectTag tag)
{
#ifdef BACKWARD_COMPATIBILITY
    if ( in.formatID() < 40 )
        in.readUInt16();
#endif
    
    ObjectTag g;
    size_t n = in.readUInt16();
    asSolid = Solid::toSolid(sim.readReference(in, g));
    Organizer::readOrganized(in, sim, n-1);
    assert_true( nbOrganized() > 0 );
    
    if ( !Buddy::check(solid()) )
        Buddy::connect(solid());
    Solid const* sol = solid();
    if ( sol->nbPoints() > 1 )
        asRadius = ( sol->posPoint(0) - sol->posPoint(1) ).norm();
    
#ifdef BACKWARD_COMPATIBILITY
    // usual number of fiber links:
    size_t nbf = nbOrganized();
    if ( in.formatID() > 50 )
#else
    size_t
#endif
        nbf = in.readUInt16();
    asLinks.resize(nbf);
    
    for ( size_t i = 0; i < nbf; ++i )
    {
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 47 )
        {
            asLinks[i].readOldFormat(in, sol);
            continue;
        }
#endif
        asLinks[i].read(in, asRadius);
        if ( asLinks[i].prime_ + DIM >= sol->nbPoints() )
            throw InvalidIO("out-of-range AsterLink index");
    }
    
    if ( nbf > 0 )
    {
        const size_t ref = asLinks[0].prime_;
        asRadius = ( sol->posPoint(ref) - sol->posPoint(ref) ).norm();
    }
}


//------------------------------------------------------------------------------
#pragma mark - Display

Vector Aster::posSolid1(size_t inx) const
{
    AsterLink const& link = asLinks[inx];
    
#ifdef BACKWARD_COMPATIBILITY
    if ( link.alt_ > 0 )
        return solid()->posPoint(link.prime_);
#endif
    
    return solid()->interpolatePoints(link.prime_, link.coef1_, link.rank_);
}

Vector Aster::posSolid2(size_t inx) const
{
    AsterLink const& link = asLinks[inx];
    
#ifdef BACKWARD_COMPATIBILITY
    if ( link.alt_ > 0 )
        return solid()->posPoint(link.alt_);
#endif
    
    return solid()->interpolatePoints(link.prime_, link.coef2_, DIM+1);
}


Vector Aster::posFiber2(size_t inx) const
{
    Fiber const* fib = fiber(inx);
    real len = asLinks[inx].len_;
    
    if ( fib->length() >= len )
    {
        return fib->posFrom(len, prop->focus);
    }
    else
    {
        // link the opposite end to an interpolation of the two solid-points:
        return fib->posEnd( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
    }
}


/**
 retrieve link between Solid and ends of Fiber
 this is only meaningfull if ( inx < nbFibers() )
 */
real Aster::getLink1(size_t inx, Vector& pos1, Vector& pos2) const
{
    pos1 = posSolid1(inx);
    if ( fiber(inx) )
    {
        pos2 = posFiber1(inx);
        return prop->stiffness[0];
    }
    pos2 = pos1;
    return 0;
}


/**
 retrieve link between Solid and side of Fiber
 this is only meaningfull if ( inx < nbFibers() )
 */
real Aster::getLink2(size_t inx, Vector& pos1, Vector& pos2) const
{
    Fiber const* fib = fiber(inx);
    
    if ( fib )
    {
        real len = asLinks[inx].len_;
        if ( fib->length() >= len )
        {
            pos1 = posSolid2(inx);
        }
        else
        {
            // interpolate between the two solid-points:
            real c = fib->length() / len;
            pos1 = ( 1.0 - c ) * posSolid1(inx) + c * posSolid2(inx);
        }
        pos2 = posFiber2(inx);
        return prop->stiffness[1];
    }

    pos1 = posSolid2(inx);
    pos2 = pos1;
    return 0;
}


/**
 This sets 'pos1' and 'pos2' as the ends of the link number `inx`
 or returns zero if the link does not exist
 */
bool Aster::getLink(size_t inx, Vector& pos1, Vector& pos2) const
{
    if ( inx < 2 * asLinks.size() )
    {
        if ( inx & 1 )
            getLink2(inx/2, pos1, pos2);
        else
            getLink1(inx/2, pos1, pos2);
        return 1;
    }
    return 0;
}


