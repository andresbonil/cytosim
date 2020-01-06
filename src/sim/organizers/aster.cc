// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

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
    assert_true( linked() );
    
    Simul & sim = simul();

    // nucleation:
    for ( size_t ii = 0; ii < asLinks.size(); ++ii )
    {
        if ( !fiber(ii) &&  RNG.test(prop->fiber_prob) )
        {
            Glossary opt(prop->fiber_spec);
            sim.add(makeFiber(sim, ii, prop->fiber_type, opt));
            opt.warnings(std::cerr, 1, " in aster:nucleate[1]");
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

/*
 Note on possible optimization:
 The coefficients of the interpolations to the Solid points are constant in time,
 and so we could simply set a matrix once, and keep it over time.
 Specifically, we would introduce a new matrix in Meca, `mK` and set it only once.
 We can then include these additional terms directly as we calculate forces in Meca:

     Y <- Y + ( mB + mC + mK ) * X

 */
void Aster::setInteractions(Meca & meca) const
{
    assert_true( linked() );

    Solid const* sol = solid();
    
    if ( !sol )
        return;

    for ( size_t n = 0 ; n < asLinks.size(); ++n )
    {
        Fiber * fib = fiber(n);

        if ( fib )
        {
            AsterLink const& link = asLinks[n];
            
            if ( link.ord == 0 )
                continue;
            
            unsigned off = sol->matIndex() + link.ref;
            unsigned pts[] = { off, off+1, off+2, off+3 };

#ifdef BACKWARD_COMPATIBILITY
            if ( link.alt > 0 )
            {
                meca.addLink(Mecapoint(sol, link.ref), fib->exactEnd(prop->focus), prop->stiffness[0]);
                if ( fib->length() > link.len )
                {
                    meca.addLink(Mecapoint(sol, link.alt), fib->interpolate(link.len, prop->focus), prop->stiffness[1]);
                }
                else
                {
                    FiberEnd tip = ( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
                    // link the opposite end to an interpolation of the two solid-points:
                    real c = fib->length() / link.len;
                    meca.addLink(fib->exactEnd(tip), Interpolation(sol, link.ref, link.alt, c), prop->stiffness[1]);
                }
                continue;
            }
#endif
            if ( link.ord == 1 )
                meca.addLink(fib->exactEnd(prop->focus), Mecapoint(sol, link.ref), prop->stiffness[0]);
            else
                meca.ADDLINK(fib->exactEnd(prop->focus), pts, link.coef1, prop->stiffness[0]);
            
            
            // make second type of link:
            real len = link.len;
            
            if ( fib->length() >= len )
            {
                if ( len > 0 )
                    meca.ADDLINK(fib->interpolate(len, prop->focus), pts, link.coef2, prop->stiffness[1]);
                else
                    meca.ADDLINK(fib->exactEnd(prop->focus), pts, link.coef2, prop->stiffness[1]);
            }
            else
            {
                // link the opposite fiber end to a new interpolation:
                FiberEnd end = ( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
                real c = fib->length() / len;
                real u = 1.0 - c;
                real coef[4];
                for ( int d = 0; d < 4; ++d )
                    coef[d] = u * link.coef1[d] + c * link.coef2[d];
                meca.ADDLINK(fib->exactEnd(end), pts, coef, prop->stiffness[1]);
            }
        }
    }
}


Aster::~Aster()
{
    //Cytosim::log("destroying %c%lu\n", TAG, identity());
    prop = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -

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
 
 */
ObjectList Aster::build(Glossary& opt, Simul& sim)
{
    assert_true(prop);
    assert_true(nbOrganized()==0);
    
    // get number of fibers:
    unsigned nbf = 0;
    std::string type, spec;
    
    opt.set(asRadius, "radius");
    opt.set(nbf,      "fibers");
    opt.set(type,     "fibers", 1);
    opt.set(spec,     "fibers", 2);
    
    Glossary fiber_opt(spec);
    
#ifdef BACKWARD_COMPATIBILITY
    if ( type.empty() && opt.set(nbf, "nb_fibers") )
    {
        type = prop->fiber_type;
        fiber_opt = opt;
        spec = "unknown";
    }
#endif
    nbOrganized(1+nbf);

    if ( asRadius <= 0 )
        throw InvalidParameter("aster:radius must be specified and > 0");

    unsigned origin = 0;
    ObjectList res = makeSolid(sim, opt, origin);
    
    if ( !solid() )
        throw InvalidParameter("could not make aster:solid");
    
    // fibers anchor points can be specified directly:
    unsigned cnt = 0;
    Vector pos1, pos2;
    std::string var = "fiber1";
    while ( opt.set(pos1, var) && opt.set(pos2, var, 1) )
    {
        //std::clog << "direct fiber anchor " << pos1 << " and " << pos2 << "\n";
        placeAnchor(pos1, pos2, origin);
        ++cnt;
        var = "fiber" + std::to_string(cnt+1);
    }

    if ( cnt < nbf )
        placeAnchors(opt, origin, nbf-cnt);
    
    //solid()->write(std::clog);
    
    if ( !spec.empty() )
    {
        for ( size_t n = 0; n < asLinks.size(); ++n )
            res.append(makeFiber(sim, n, type, fiber_opt));
    }
    
    return res;
}


ObjectList Aster::makeFiber(Simul& sim, size_t inx, std::string const& fiber_type, Glossary& opt)
{
    ObjectList objs = sim.fibers.newObjects(fiber_type, opt);
    
    if ( objs.empty() )
        throw InvalidParameter("could not create aster:fiber");

    Fiber * fib = Fiber::toFiber(objs[0]);

    if ( !fib )
        throw InvalidParameter("unexpected object returned by fibers.newObjects()");

    grasp(fib, inx+1);

    Vector pos = posLink1(inx);
    Vector dir = posLink2(inx) - pos;
    real n = dir.normSqr();
    
    if ( n > REAL_EPSILON )
    {
        if ( prop->focus == PLUS_END )
            dir /= -sqrt(n);
        else
            dir /= sqrt(n);
    }
    else
        dir = Vector::randU();
    
    //std::clog << "new aster:fiber " << pos << " and " << dir << "\n";
    ObjectSet::rotateObjects(objs, Rotation::rotationToVector(dir));
    ObjectSet::translateObjects(objs, pos - fib->posEnd(prop->focus));
    
    return objs;
}


/**
 Anchor a Fiber between positions A and B, specified in a local reference frame
 associated with the Aster. Dimensions will be scaled by 'asRadius' eventually.
 */
void Aster::placeAnchor(Vector const& A, Vector const& B, unsigned ref)
{
    AsterLink & link = asLinks.new_val();
    //std::clog << "Aster::placeAnchor(" << asLinks.size() << ")\n";
    link.set(A, B);
    link.len *= asRadius;
    link.ref = ref;
    //link.print(std::clog);
}


ObjectList Aster::makeSolid(Simul& sim, Glossary& opt, unsigned& origin)
{
    ObjectList res(2);
    Solid * sol = nullptr;
    
    // find the Solid specified:
    std::string spec;
    if ( opt.set(spec, "solid") )
    {
        SolidProp * p = static_cast<SolidProp*>(sim.findProperty("solid", spec));
        
        if ( p )
        {
            sol = new Solid(p);
            res.push_back(sol);
            res.append(sol->build(opt, sim));
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
        }
    }
    
    if ( ! sol )
        throw InvalidParameter("aster:solid must be specified");

    // check that there is at least one point:
    if ( sol->dragCoefficient() < REAL_EPSILON )
#ifdef BACKWARD_COMPATIBILITY
        sol->addSphere(Vector(0,0,0), asRadius);
#else
        throw InvalidParameter("Aster's drag coefficient is null: please specify 'point1=center, RADIUS'");
#endif
    
    // add local coordinate system around the last point:
    origin = sol->addTriad(asRadius);
    sol->fixShape();
    
    //std::cerr << *sol << '\n';
    grasp(sol, 0);
    return res;
}


/**
 One can specify the `radius` of the aster, and `nb_fibers`.
 
 The aster 'type' can be:
 - `astral` fiberd are anchored at random positions near the center, pointing outward
 - `radial` fibers are anchored always at the same distance from the center, pointing radially
 - `regular` fibers are anchored regularly over the surface and point radially
 - `angular` where all fibers are restricted within an specified solid angle,
 .
 */
void Aster::placeAnchors(Glossary & opt, unsigned origin, unsigned nbf)
{
    real dis = 0;
    if ( opt.set(dis, "radius", 1) && dis > asRadius )
        throw InvalidParameter("aster:radius[1] must be smaller than aster:radius[0]");

    const real alpha = dis / asRadius;
    
    unsigned type = 0;
    opt.set(type, "type", {{"radial", 0}, {"astral", 1}, {"regular", 2}, {"angular", 3}, {"disc", 4}});
    
    if ( type == 4 )
    {
        // This is a special case for Yeast spindles
        // use a separation of 25 nm by default, corresponding to Microtubules
        real sep = 0.025;
        opt.set(sep, "seed_diameter");
        std::vector<Vector2> pts(nbf);
        size_t ouf = 0;
        size_t cnt = 0;
        do {
            cnt = tossPointsDisc(pts, sep/asRadius, 1024);
        } while ( cnt < nbf && ++ouf < 1024 );
        if ( cnt < nbf )
        {
            std::clog << "warning: aster could only fit " << cnt << " seeds ";
            std::clog << "with aster:seed_diameter = " << sep << '\n';
        }
        //std::clog << "toss(" << nbf << ") placed " << cnt << "\n";
        real d = dis * 0.5;
        for ( size_t n = 0; n < cnt; ++n )
        {
            // orient anchors by default along the X-axis:
            real y = pts[n].YY;
#if ( DIM == 2 )
            placeAnchor(Vector2(-d,y), Vector2(+d,y), origin);
#elif ( DIM == 3 )
            real x = pts[n].XX;
            placeAnchor(Vector3(-d,x,y), Vector3(+d,x,y), origin);
#endif
        }
    }
    else if ( type == 3 )
    {
        /*
         For type 'angular' all fibers are restricted within an specified solid angle,
         and their orientation is radial
         by GAELLE LETORT, 14.03.2017
         */
        real angle = M_PI;
        opt.set(angle, "angle") || opt.set(angle, "aster_angle");
        if ( angle < REAL_EPSILON )
            throw InvalidParameter("aster::angle must be > 0");
#if ( DIM == 1 )
        // No effect of angle in 1D, same as default
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector D(n%2?1:-1);
            placeAnchor(Vector(0.0), D, origin);
        }
#elif ( DIM == 2 )
        real delta = 2 * angle / real(nbf);
        // points are evenly distributed from -aster_angle to aster_angle
        real ang = -angle;
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector P(cos(ang), sin(ang));
            placeAnchor(alpha*P, P, origin);
            ang += delta;
        }
#else
        real cap = 1.0 - cos(angle);
        real sep, sep0 = sqrt( 2 * M_PI * cap / nbf );
        std::vector<Vector> pts(nbf, Vector(0,0,0));
        size_t ouf = 0;
        size_t cnt = 0;
        do {
            // we decrease gradually the separation, to reach a good solution...
            sep = 512 * sep0 / real(ouf+512);
            cnt = tossPointsCap(pts, cap, sep, 1024);
            //std::clog << "toss(" << nbf << ") placed " << cnt << " with sep = " << sep << "\n";
        } while ( cnt < nbf && ++ouf < 1024 );
        if ( cnt < nbf )
            std::clog << "warning: aster could only fit " << cnt << " seeds\n";
        //std::clog << "toss(" << nbf << ") placed " << cnt << " with sep = " << sep << "\n";
        for ( size_t n = 0; n < cnt; ++n )
            placeAnchor(alpha*pts[n], pts[n], origin);
#endif
    }
    else if ( type == 2 )
    {
        /*
         For type 'regular' we put fibers regularly on the surface,
         */
#if ( DIM == 1 )
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector D(n%2?1:-1);
            placeAnchor(Vector(0.0), D, origin);
        }
#elif ( DIM == 2 )
        real ang = 0, delta = 2 * M_PI / real(nbf);
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector P(cos(ang), sin(ang));
            placeAnchor(alpha*P, P, origin);
            ang += delta;
        }
#else
        //we use PointsOnSphere to distribute points 'equally' on the sphere
        PointsOnSphere sphere(nbf);
        Vector P;
        for ( unsigned n = 0; n < nbf; ++n )
        {
            sphere.copyPoint(P, n);
            placeAnchor(alpha*P, P, origin);
        }
#endif
    }
    else if ( type == 1 )
    {
        /*
         For type 'astral' we put fibers randomly on the surface,
         with a constrain based on the scalar product: position*direction > 0
         */
        if ( dis <= 0 )
            throw InvalidParameter("aster:radius[1] must be specified and >= 0");
        for ( unsigned n = 0; n < nbf; ++n )
        {
            Vector P = Vector::randB();
            Vector D = Vector::randU();
            while ( dot(D, P) < 0 )
                D = Vector::randU();
            placeAnchor(P-alpha*D, P, origin);
        }
    }
    else if ( type == 0 )
    {
        /*
         For type 'radial' we put fibers randomly on the surface, and set their
         direction as purely radial. We require a separation of 25 nm by default,
         corresponding to Microtubule's size.
         */
        real sep = 0.025;
        opt.set(sep, "seed_diameter");
        size_t ouf = 0;
        size_t cnt = 0;
        std::vector<Vector> pts(nbf, Vector(0,0,0));
        do {
            cnt = tossPointsSphere(pts, sep/asRadius, 1024);
        } while ( cnt < nbf && ++ouf < 1024 );
        if ( cnt < nbf )
        {
            std::clog << "warning: aster could only fit " << cnt << " seeds ";
            std::clog << "with aster:seed_diameter = " << sep << std::endl;
        }
        //std::clog << "toss(" << nbf << ") placed " << cnt << "\n";
        for ( size_t n = 0; n < cnt; ++n )
            placeAnchor(alpha*pts[n], pts[n], origin);
    }
    else
        throw InvalidParameter("unknown aster:type");
}


//------------------------------------------------------------------------------
#pragma mark -

void Aster::write(Outputter& out) const
{
    Organizer::write(out);
    
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
    
    Organizer::read(in, sim, tag);
    
    assert_true( nbOrganized() > 0 );
    assert_true( organized(0)->tag() == Solid::TAG );
    
    Solid * sol = solid();
    if ( sol->nbPoints() > 1 )
        asRadius = ( sol->posPoint(0) - sol->posPoint(1) ).norm();
    
#ifdef BACKWARD_COMPATIBILITY
    // usual number of fiber links:
    unsigned nbf = nbOrganized() - 1;
    if ( in.formatID() > 50 )
#else
    unsigned
#endif
        nbf = in.readUInt16();
    asLinks.resize(nbf);
    
    for ( size_t i = 0; i < nbf; ++i )
    {
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 47 )
        {
            asLinks[i].reset();
            unsigned a = in.readUInt16();
            unsigned b = in.readUInt16();
            if ( b >= sol->nbPoints() )
                throw InvalidIO("invalid AsterLink index");
            asLinks[i].ref = a;
            asLinks[i].coef1[0] = 1.0;
            asLinks[i].alt = b;
            asLinks[i].len = (sol->posPoint(a)-sol->posPoint(b)).norm();
            continue;
        }
#endif
        asLinks[i].read(in);
        asLinks[i].len *= asRadius;
        if ( asLinks[i].ref + asLinks[i].ord >= sol->nbPoints() )
            throw InvalidIO("invalid AsterLink index");
    }
    
    if ( nbf > 0 )
    {
        unsigned ref = asLinks[0].ref;
        asRadius = ( sol->posPoint(ref) - sol->posPoint(ref) ).norm();
    }
}


//------------------------------------------------------------------------------
#pragma mark -

Vector Aster::posLink1(size_t inx) const
{
    Solid const* sol = solid();
    real const* coef = asLinks[inx].coef1;
    const unsigned ref = asLinks[inx].ref;
    
#ifdef BACKWARD_COMPATIBILITY
    if ( asLinks[inx].alt > 0 )
        return sol->posPoint(ref);
#endif

    int top = std::min(DIM+1u, sol->nbPoints());
    Vector res = coef[0] * sol->posPoint(ref);
    for ( int i = 1; i < top; ++i )
        res += coef[i] * sol->posPoint(i+ref);
    
    return res;
}

Vector Aster::posLink2(size_t inx) const
{
    Solid const* sol = solid();
    real const* coef = asLinks[inx].coef2;
    const unsigned ref = asLinks[inx].ref;
    
#ifdef BACKWARD_COMPATIBILITY
    if ( asLinks[inx].alt > 0 )
        return sol->posPoint(asLinks[inx].alt);
#endif

    int top = std::min(DIM+1u, sol->nbPoints());
    Vector res = coef[0] * sol->posPoint(ref);
    for ( int i = 1; i < top; ++i )
        res += coef[i] * sol->posPoint(i+ref);
    
    return res;
}

Vector Aster::posFiber2(size_t inx) const
{
    Fiber const* fib = fiber(inx);
    real len = asLinks[inx].len;
    
    if ( fib->length() >= len )
    {
        if ( len > 0 )
            return fib->pos(len, prop->focus);
        else
            return fib->posEnd(prop->focus);
    }
    else
    {
        // link the opposite end to an interpolation of the two solid-points:
        return fib->posEnd( prop->focus == PLUS_END ? MINUS_END : PLUS_END );
    }
}

/**
 This sets 'pos1' and 'pos2' as the ends of the link number `inx`
 or returns zero if the link does not exist
 */
bool Aster::getLink(size_t inx, Vector& pos1, Vector& pos2) const
{
    size_t n = inx / 2;
    
    if ( n < asLinks.size() && asLinks[n].ord > 0 )
    {
        Fiber const* fib = fiber(n);
        
        if ( inx & 1 )
        {
            pos1 = posLink1(n);
            if ( fib )
                pos2 = fib->posEnd(prop->focus);
            else
                pos2 = pos1;
        }
        else
        {
            if ( fib )
            {
                real len = asLinks[n].len;
                if ( fib->length() >= len )
                {
                    pos1 = posLink2(n);
                }
                else
                {
                    // interpolate between the two solid-points:
                    real c = fib->length() / len;
                    pos1 = ( 1.0 - c ) * posLink1(n) + c * posLink2(n);
                }
                pos2 = posFiber2(n);
            }
            else
            {
                pos1 = posLink2(n);
                pos2 = pos1;
            }
        }
        return true;
    }
    return false;
}


