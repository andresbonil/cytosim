// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "glossary.h"
#include "mecapoint.h"
#include "organizer_set.h"
#include "nucleus.h"
#include "bundle.h"
#include "aster.h"
#include "fake.h"
#include "solid.h"
#include "simul.h"

//------------------------------------------------------------------------------

void OrganizerSet::step()
{
    Organizer * obj = first(), * nxt;
    while ( obj )
    {
        nxt = obj->next();
        obj->step();
        obj = nxt;
    }
}

//------------------------------------------------------------------------------

Property* OrganizerSet::newProperty(const std::string& cat, const std::string& nom, Glossary&) const
{
    if ( cat == "aster" )   return new AsterProp(nom);
    if ( cat == "bundle" )  return new BundleProp(nom);
    if ( cat == "nucleus" ) return new NucleusProp(nom);
    if ( cat == "fake" )    return new FakeProp(nom);
    return nullptr;
}


Object * OrganizerSet::newObjectT(const ObjectTag tag, unsigned num)
{
    if ( tag == Aster::TAG )
    {
        AsterProp * p = simul.findProperty<AsterProp>("aster", num);
        return new Aster(p);
    }
    
    if ( tag == Bundle::TAG )
    {
        BundleProp * p = simul.findProperty<BundleProp>("bundle", num);
        return new Bundle(p);
    }
    
    if ( tag == Nucleus::TAG )
    {
        NucleusProp * p = simul.findProperty<NucleusProp>("nucleus", num);
        return new Nucleus(p);
    }
    
    if ( tag == Fake::TAG )
    {
        FakeProp * p = simul.findProperty<FakeProp>("fake", num);
        return new Fake(p);
    }
    
    return nullptr;
}


ObjectList OrganizerSet::newObjects(const std::string& nom, Glossary& opt)
{
    Organizer * obj = nullptr;
    Property * p = simul.properties.find_or_die(nom);
    
    if ( p->category() == "aster" )
        obj = new Aster(static_cast<AsterProp*>(p));
    else if ( p->category() == "bundle" )
        obj = new Bundle(static_cast<BundleProp*>(p));
    else if ( p->category() == "nucleus" )
        obj = new Nucleus(static_cast<NucleusProp*>(p));
    else if ( p->category() == "fake" )
        obj = new Fake(static_cast<FakeProp*>(p));

    ObjectList res;
    if ( obj )
    {
        res = obj->build(opt, simul);
        res.push_back(obj);
    }
    
    return res;
}

//------------------------------------------------------------------------------

void OrganizerSet::add(Object * obj)
{
    ObjectSet::add(obj);
    // we also link all dependent objects:
    static_cast<Organizer*>(obj)->addOrganized(simul);
}


Aster * OrganizerSet::findAster(const ObjectID n) const
{
    return Aster::toAster(findID(n));
}


Organizer * OrganizerSet::findOrganizer(const Mecable * m) const
{
    for ( Organizer * o=first(); o; o=o->next() )
        for ( unsigned i = 0; i < o->nbOrganized(); ++i )
            if ( m == o->organized(i) )
                return o;

    return nullptr;
}


//------------------------------------------------------------------------------
void OrganizerSet::foldPosition(Modulo const* s) const
{
    for ( Organizer * o=first(); o; o=o->next() )
        o->foldPosition(s);
}


void OrganizerSet::report(std::ostream& os) const
{
    if ( size() > 0 )
    {
        unsigned total = 0;
        os << title() << "\n";
        for ( Property * i : simul.properties.find_all("aster", "bundle") )
        {
            unsigned cnt = count(match_property, i);
            os << std::setw(10) << cnt << " " << i->name() << '\n';
            ++total;
        }
        for ( Property * i : simul.properties.find_all("nucleus", "fake") )
        {
            unsigned cnt = count(match_property, i);
            os << std::setw(10) << cnt << " " << i->name() << '\n';
            ++total;
        }
        if ( total > 1 )
            os << std::setw(10) << size() << " total\n";
    }
}

