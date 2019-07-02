// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "sim.h"
#include "simul.h"
#include "exceptions.h"
#include "hand_prop.h"
#include "simul_prop.h"
#include "backtrace.h"
#include "modulo.h"
#include "tictoc.h"

extern Modulo const* modulo;

#include "simul_step.cc"
#include "simul_file.cc"
#include "simul_custom.cc"
#include "simul_report.cc"
#include "simul_solve.cc"

#include "nucleus.h"
#include "aster.h"
#include "fake.h"
#include "wrist.h"
#include "bundle.h"
#include "space_strip.h"
#include "space_periodic.h"
#include "space_cylinderP.h"
#include "fiber.h"
#include "event.h"

#include <csignal>

//---------------------------  global variables/functions ---------------------

void out_of_memory_handler()
{
    write(STDERR_FILENO, "\n* * * * *\n", 11);
    write(STDERR_FILENO, "Cytosim: memory allocation failed", 33);
    write(STDERR_FILENO, "\n* * * * *\n", 11);
    print_backtrace();
    _exit(1);
}

void termination_handler()
{
    write(STDERR_FILENO, "\n* * * * *\n", 11);
    write(STDERR_FILENO, "Cytosim: uncaught exception", 27);
    write(STDERR_FILENO, "\n* * * * *\n", 11);
    print_backtrace();
    abort();
}

void signal_handler(int sig)
{
    write(STDERR_FILENO, "\n* * * * *\n", 11);
    psignal(sig, "Cytosim");
    write(STDERR_FILENO, "* * * * *\n", 10);
    print_backtrace();
    _exit(sig);
}

//------------------------------------------------------------------------------
#pragma mark -

Simul::Simul()
: prop(nullptr), spaces(*this), fields(*this),
fibers(*this), spheres(*this), beads(*this), solids(*this),
singles(*this), couples(*this), organizers(*this), events(*this)
{
    sReady        = false;
    precondCPU[0] = 0;
    precondCPU[1] = 0;
    precondCPU[2] = 0;
    precondCPU[3] = 0;
    precondMethod = 1;
    precondCounter = 0;
    
    prop = new SimulProp("undefined");
}

Simul::~Simul()
{
    erase();
    
    if ( prop )
        delete(prop);
}

//------------------------------------------------------------------------------
#pragma mark -

Space const* Simul::findSpace(std::string const& str) const
{
    if ( str == "first" )
        return static_cast<Space*>(spaces.inventory.first());

    if ( str == "last" )
        return static_cast<Space*>(spaces.inventory.last());
    
    Property * sp = properties.find("space", str);
    
    if ( sp )
        return spaces.findObject(sp);
    else
        return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This will initialize the simulation by registering callbacks.
 You should still call Simul::prepare() before calling step()
 */
void Simul::initialize(Glossary & glos)
{
    // Register a function to be called if operator new fails:
    std::set_new_handler(out_of_memory_handler);
    
    // Register a function to be called upon abortion:
    std::set_terminate(termination_handler);
    
    // Register a function to be called for Floating point exceptions:
    if ( signal(SIGFPE, signal_handler) == SIG_ERR )
        std::cerr << "Could not register SIGFPE handler\n";
    
    // read parameters, and complete
    prop->read(glos);
}

//------------------------------------------------------------------------------

real Simul::time() const
{
    return prop->time;
}

void Simul::erase()
{
    //std::cerr << "Simul::erase()\n";
    relax();
    organizers.erase();
    fibers.erase();
    singles.erase();
    couples.erase();
    spheres.erase();
    beads.erase();
    solids.erase();
    fields.erase();
    spaces.erase();
    events.erase();
    
    // destroy all properties, except the SimulProp:
    properties.erase();
 
    prop->time = 0;
    modulo     = nullptr;
}


size_t Simul::nbObjects() const
{
    return  (  organizers.size()
             + singles.size()
             + couples.size()
             + fibers.size()
             + beads.size()
             + solids.size()
             + spheres.size()
             + spaces.size()
             + fields.size() );
}


void Simul::foldPosition() const
{
    if ( modulo )
    {
        fibers.foldPosition(modulo);
        beads.foldPosition(modulo);
        solids.foldPosition(modulo);
        spheres.foldPosition(modulo);
        singles.foldPosition(modulo);
        couples.foldPosition(modulo);
        organizers.foldPosition(modulo);
    }
}


//------------------------------------------------------------------------------
#pragma mark -

/**
 Convert Object pointer to Mecable* if possible
 */
Mecable* Simul::toMecable(Object * obj)
{
    if ( obj )
    switch( obj->tag() )
    {
        case  Fiber::TAG:  return static_cast<Mecable*>(obj);
        case   Bead::TAG:  return static_cast<Mecable*>(obj);
        case  Solid::TAG:  return static_cast<Mecable*>(obj);
        case Sphere::TAG:  return static_cast<Mecable*>(obj);
    }
    return nullptr;
}

/**
 Find an object from a Human-friendly representation, such as
 fiber1
 single1
 */
Mecable * Simul::findMecable(const std::string& arg) const
{
    Object  * obj = fibers.findObject(arg);
    if (!obj) obj = solids.findObject(arg);
    if (!obj) obj = spheres.findObject(arg);
    if (!obj) obj = beads.findObject(arg);
    return static_cast<Mecable*>(obj);
}


void Simul::add(Object * w)
{
    //std::clog << " Simul::add(" << w->reference() << ")" << std::endl;
    assert_true(w);
    ObjectSet * set = findSetT(w->tag());
    set->add(w);
}


void Simul::add(ObjectList const& objs)
{
    //std::clog << " Simul::add("<< objs.size() <<" objects):" << std::endl;
    for ( Object * obj : objs )
        if ( obj )
            add(obj);
}


void Simul::remove(Object * w)
{
    assert_true( w->objset() );
    w->objset()->remove(w);
}


void Simul::remove(ObjectList const& objs)
{
    //std::clog << " Simul::remove("<< objs.size() <<" objects):" << std::endl;
    for ( Object * obj : objs )
        if ( obj )
            remove(obj);
}


void Simul::erase(Object * w)
{
    //std::clog << "Simul::erase " << w->reference() << std::endl;
    remove(w);
    delete(w);
}


void Simul::erase(ObjectList const& objs)
{
    //std::clog << " Simul::erase(" << objs.size() << " objects):\n";
    for ( Object * obj : objs )
        if ( obj )
        {
            //std::clog << " Simul::erase(" << obj << ")\n";
            remove(obj);
            delete(obj);
        }
}


void Simul::mark(ObjectList const& objs, ObjectMark mrk)
{
    //std::clog << " Simul::erase("<< objs.size() <<" objects):" << std::endl;
    for ( Object * i : objs )
        i->mark(mrk);
}

//------------------------------------------------------------------------------
#pragma mark -

/**
 This is used primarily to parse the configuration file,
 using full class name
 */
ObjectSet * Simul::findSet(const std::string& cat)
{
    //std::clog << "findSet("<<kind<<")"<<std::endl;
    if ( cat == "space" )        return &spaces;
    if ( cat == "field" )        return &fields;
    if ( cat == "fiber" )        return &fibers;
    if ( cat == "bead" )         return &beads;
    if ( cat == "solid" )        return &solids;
    if ( cat == "sphere" )       return &spheres;
    if ( cat == "single" )       return &singles;
    if ( cat == "couple" )       return &couples;
    if ( cat == "organizer" )    return &organizers;
    if ( cat == "aster" )        return &organizers;
    if ( cat == "bundle" )       return &organizers;
    if ( cat == "nucleus" )      return &organizers;
    if ( cat == "fake" )         return &organizers;
    if ( cat == "event" )        return &events;
    return nullptr;
}


/**
 This is used primarily to read the binary trajectory file,
 using a single character to refer to each class in Cytosim
 */
ObjectSet * Simul::findSetT(const ObjectTag tag)
{
    switch( tag )
    {
        case        Couple::TAG:    return &couples;
        case        Single::TAG:    return &singles;
        case         Wrist::TAG:    return &singles;
        case         Fiber::TAG:    return &fibers;
        case Fiber::TAG_DYNAMIC:    return &fibers;
        case Fiber::TAG_LATTICE:    return &fibers;
        case          Bead::TAG:    return &beads;
        case         Solid::TAG:    return &solids;
        case        Sphere::TAG:    return &spheres;
        case         Field::TAG:    return &fields;
        case         Space::TAG:    return &spaces;
        case       Nucleus::TAG:    return &organizers;
        case        Bundle::TAG:    return &organizers;
        case         Aster::TAG:    return &organizers;
        case          Fake::TAG:    return &organizers;
        case         Event::TAG:    return &events;
        case        Object::TAG:    return nullptr;
#ifdef BACKWARD_COMPATIBILITY
        case 'm':                   return &fibers;
#endif
    }
    return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark -

bool Simul::isPropertyClass(const std::string& name) const
{
    if ( name == "simul" )
        return true;
    
    if ( name == "hand" )
        return true;
    
    return const_cast<Simul*>(this)->findSet(name);
}


Property* Simul::findProperty(const std::string& cat, const std::string& nom) const
{
    if ( cat == "simul" && nom == prop->name() )
        return prop;

    return properties.find(cat, nom);
}


Property* Simul::findProperty(const std::string& nom) const
{
    if ( nom == prop->name() )
        return prop;
    
    return properties.find(nom);
}


PropertyList Simul::findAllProperties(const std::string& cat) const
{
    if ( cat == "simul" )
    {
        PropertyList list;
        list.push_back(prop);
        return list;
    }
    
    return properties.find_all(cat);
}


/**
 @defgroup ObjectGroup List of objects
 
 The command `set simul` will define the global parameters.
 The `simul` is automatically created, and you cannot use 'new simul'.

 Objects       | Base class    | Parameters    |
 --------------|---------------|----------------
 `simul`       |  Simul        | @ref SimulPar  
 
 
 These objects cannot move:
 
 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `space`       |  Space        | @ref SpacePar    | @ref SpaceGroup
 `field`       |  Field        | @ref FieldPar    | -                 
 `event`       |  Event        | @ref EventPar    | -
 
 
 `Mecables` can move or deform, and come in 4 basic forms:
 
 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `fiber`       |  Fiber        | @ref FiberPar    | @ref FiberGroup
 `bead`        |  Bead         | @ref SolidPar    | -
 `solid`       |  Solid        | @ref SolidPar    | -
 `sphere`      |  Sphere       | @ref SpherePar   | -

 
 A `Hand` is an object that can bind to fiber, but it can only be used
 as a sub-part of `Single` or `Couple`.

 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `hand`        |  Hand         | @ref HandPar     | @ref HandGroup
 
 
 `Single` and `Couple` contain one or two `Hand` respectively:

 Class Name    | Base class    | Parameters       | Specialization   |
 --------------|---------------|------------------|-------------------
 `single`      |  Single       | @ref SinglePar   | @ref SingleGroup
 `couple`      |  Couple       | @ref CouplePar   | @ref CoupleGroup
 
 
 
 The `Organizers` describe composite objects build from multiple Mecables:
 
 Organizers    | Base class    | Parameters      |
 --------------|---------------|------------------
 `aster`       |  Aster        | @ref AsterPar    
 `bundle`      |  Bundle       | @ref BundlePar   
 `nucleus`     |  Nucleus      | @ref NucleusPar  
 `fake`        |  Fake         | @ref FakePar     
 .
 
 */
Property* Simul::newProperty(const std::string& cat, const std::string& nom, Glossary& glos)
{
    if ( cat == "simul" )
    {
        /* There can only be one SimulProp and it is already created */
        assert_true(prop);
        prop->rename(nom);
        //std::clog << "Simul is named `" << nom << "'\n";
        return prop;
    }
    
    if ( isPropertyClass(nom) )
        throw InvalidSyntax("`"+nom+"' is a reserved keyword");
    
    Property * p = findProperty(nom);
    
    if ( p )
        throw InvalidSyntax("property `"+nom+"' is already defined");
    
    if ( cat == "hand" )
    {
        p = HandProp::newProperty(nom, glos);
        properties.deposit(p);
    }
    else
    {
        ObjectSet * set = findSet(cat);
        
        if ( !set )
            throw InvalidSyntax("unknown class `"+cat+"'");
        
        p = set->newProperty(cat, nom, glos);
        properties.deposit(p);
    }
    
    return p;
}
