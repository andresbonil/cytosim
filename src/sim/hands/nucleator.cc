// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "nucleator.h"
#include "nucleator_prop.h"
#include "glossary.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "fiber_prop.h"
#include "fiber_set.h"
#include "hand_monitor.h"
#include "simul.h"


//------------------------------------------------------------------------------

Nucleator::Nucleator(NucleatorProp const* p, HandMonitor* h)
: Hand(p,h), prop(p)
{
    gspTime = RNG.exponential();
}

//------------------------------------------------------------------------------

void Nucleator::makeFiber(Simul& sim, Vector pos, std::string const& fiber_type, Glossary& opt)
{    
    ObjectList objs = sim.fibers.newObjects(fiber_type, opt);
    if ( objs.empty() )
        return;

    Fiber * fib = Fiber::toFiber(objs[0]);
    
    // register the new objects:
    sim.add(objs);
    
    // indicate the origin of nucleation:
    ObjectMark mk = 0;
    if ( opt.set(mk, "mark") )
        Simul::mark(objs, mk);
    else
        Simul::mark(objs, haMonitor->nucleatorID());

    // the Fiber will be oriented depending on specificity:
    Rotation rot;
    
    real ang = 0;
    if ( opt.set(ang, "nucleation_angle") )
    {
        Vector dir = haMonitor->otherDirection(this);
        rot = Rotation::rotationToVector(dir);
#if ( DIM == 2 )
        real c = cos(ang), s = RNG.sflip() * sin(ang);
        rot = rot * Rotation::rotation(c, s);
#elif ( DIM == 3 )
        rot = rot * Rotation::rotationAroundX(RNG.sreal()*M_PI) * Rotation::rotationAroundZ(ang);
#endif
    }
    else switch( prop->specificity )
    {            
        case NucleatorProp::NUCLEATE_PARALLEL:
        {
            Vector dir = haMonitor->otherDirection(this);
            rot = Rotation::randomRotationToVector(dir);
        }
        break;
        
        case NucleatorProp::NUCLEATE_ANTIPARALLEL:
        {
            Vector dir = -haMonitor->otherDirection(this);
            rot = Rotation::randomRotationToVector(dir);
        }
        break;
        
        case NucleatorProp::NUCLEATE_PARALLEL_IF:
        {
            Hand * ha = haMonitor->otherHand(this);
            if ( ha && ha->attached() )
            {
                rot = Rotation::randomRotationToVector(ha->dirFiber());
                // remove key to avoid warning:
                opt.clear("orientation");
            }
            else
            {
                fib->mark(0);
                std::string str;
                if ( opt.set(str, "orientation") )
                {
                    std::istringstream iss(str);
                    rot = Movable::readRotation(iss, pos, fib->prop->confine_space_ptr);
                }
                else {
                    rot = Rotation::randomRotation();
                }
            }
        }
        break;
        
        case NucleatorProp::NUCLEATE_ORIENTATED:
        {
            std::string str;
            if ( opt.set(str, "orientation") )
            {
                std::istringstream iss(str);
                rot = Movable::readRotation(iss, pos, fib->prop->confine_space_ptr);
            }
            else {
                rot = Rotation::randomRotation();
            }
        }
        break;

        default:
            throw InvalidParameter("unknown nucleator:specificity");
    }
    
    ObjectSet::rotateObjects(objs, rot);
    
    // shift position by the length of the interaction:
    if ( haMonitor->interactionLength() > 0 )
    {
        Vector dir = haMonitor->otherDirection(this);
        pos += dir.randOrthoU(haMonitor->interactionLength());
    }

    /*
     We translate Fiber to match the Nucleator's position,
     and if prop->hold_end, the Hand is attached to the new fiber
     */
    if ( prop->hold_end == MINUS_END )
    {
        attachEnd(fib, MINUS_END);
        ObjectSet::translateObjects(objs, pos-fib->posEndM());
    }
    else if ( prop->hold_end == PLUS_END )
    {
        attachEnd(fib, PLUS_END);
        ObjectSet::translateObjects(objs, pos-fib->posEndP());
    }
    else
        ObjectSet::translateObjects(objs, pos-fib->position());

    //std::clog << "nucleated fiber in direction " << fib->dirEndM() << "\n";

    // report unused values:
    opt.warnings(std::cerr, 1, " in nucleator's spec");
}


//------------------------------------------------------------------------------
/**
 Does not attach nearby Fiber, but can nucleate
 */
void Nucleator::stepUnattached(const FiberGrid&, Vector const& pos)
{
    assert_false( attached() );
    
    gspTime -= prop->rate_dt;
    
    if ( gspTime < 0 )
    {
        gspTime = RNG.exponential();
        try {
            Glossary opt(prop->fiber_spec);
            makeFiber(*haMonitor->simul(), pos, prop->fiber_type, opt);
        }
        catch( Exception & e )
        {
            e << "\nException occurred while executing nucleator:code";
            throw;
        }
    }
}


void Nucleator::stepUnloaded()
{
    assert_true( attached() );
    
    if ( testDetachment() )
        return;
    
    /// OPTION 1: delete entire fiber
    if ( prop->addictive == 2 )
    {
        delete(fiber());
        return;
    }
    
    // may track the end of the Fiber:
    if ( prop->track_end == MINUS_END )
        relocateM();
    else if ( prop->track_end == PLUS_END )
        relocateP();
}


void Nucleator::stepLoaded(Vector const& force, real force_norm)
{
    assert_true( attached() );
    assert_true( nextDetach >= 0 );
    
    if ( testKramersDetachment(force_norm) )
        return;
    
    // may track the end of the Fiber:
    if ( prop->track_end == MINUS_END )
        relocateM();
    else if ( prop->track_end == PLUS_END )
        relocateP();
}


//------------------------------------------------------------------------------
/**
 If prop->addictive, this gives a poisonous goodbye-kiss to the fiber
 */
void Nucleator::detach()
{
    if ( prop->addictive )
        fiber()->setDynamicState(nearestEnd(), STATE_RED);
    
    Hand::detach();
}

