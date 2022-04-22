#ifndef REPORT_PYTHON_H
#define REPORT_PYTHON_H

#include <fstream>
#include <sstream>

#include "stream_func.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "splash.h"
#include "parser.h"
#include "simul.h"
#include "simul_prop.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "objecter_python.h"
#include "objecter_python.h"
#include "property_modules.h"
#include "fiber_modules.h"
#include "solid_modules.h"
#include "space_modules.h"
#include "single_modules.h"
#include "couple_modules.h"
#include "hand_modules.h"
#include "python_utilities.h"
namespace py = pybind11;

class Simul;
class SimulProp;


/// ObjGroup : a vector of objects of same type having the same property
template<typename Obj, typename Prp> 
class ObjGroup : public std::vector<Obj*>{
    public:
    Prp * prop;
    ObjGroup() = default;
    ObjGroup(Prp * p) : ObjGroup() {prop = p ;};
    ~ObjGroup() = default;
};

/// ObjMap : a map <string, ObjGroup>
template<typename Obj, typename Prp> 
using ObjMap = std::map<std::string,ObjGroup<Obj,Prp>> ;

/// A time frame ; basically a wrapper around a object dictionnary
class Frame 
{
public:
        /// An Objectmap is map of  (string,objectgroup)
        ObjMap<Fiber,FiberProp> fibers;
        ObjMap<Solid,SolidProp> solids;
        ObjMap<Space,SpaceProp> spaces;
        ObjMap<Couple,CoupleProp> couples;
        ObjMap<Single,SingleProp> singles;

        // Time of the frame
        real time;
        int index;
        
        /// pointer to simul
        Simul * simul;
        
        /// The party zone
        py::dict objects;
        
        /// Default constr and destrc
        Frame() = default;
        ~Frame() = default;
};

/// Distribute the objects (pointers) in the groups and in the dict.
template<typename Obj, typename Prp, typename Set> 
void distribute_objects(Simul * sim, Frame * current, ObjMap<Obj,Prp> mappe, Set & set, std::string categ ) {
    // First we list all objects in category, and create the ObjGroups in the map
    PropertyList plist = sim->properties.find_all(categ);
    if (!plist.empty()) {
        for ( Property * i : plist )
            {
                Prp * fp = static_cast<Prp*>(i);
                mappe[fp->name()] = ObjGroup<Obj,Prp>(fp);
            }
        // Then we assign all objects to their groups
        Obj * obj = set.first();
        while (obj) {
            mappe[obj->property()->name()].push_back(obj);
            obj = obj->next();
        }
        // Then we fill the dictionnary
        for (const auto &[name, group] : mappe) {
            current->objects[py::cast(name)] = group;
        }
        
    }
}
 
/// Distribute the objects (pointers) in the groups and in the dict ; 
// special case for couple, single, where firstID needs to be used
template<typename Obj, typename Prp, typename Set> 
void distribute_objects_wID(Simul * sim, Frame * current, ObjMap<Obj,Prp> mappe, Set & set, std::string categ )
{
    // First we list all objects in category, and create the ObjGroups in the map
    PropertyList plist = sim->properties.find_all(categ);
    if (!plist.empty()) {
        for ( Property * i : plist )
            {
                Prp * fp = static_cast<Prp*>(i);
                mappe[fp->name()] = ObjGroup<Obj,Prp>(fp);
            }
        // Then we assign all objects to their groups
        // (OUTDATED) We need to add a static cast here because ...
        //      sometimes first, last comme from the base class ObjectSet, 
        //      sometimes from a derived class, e.g. FiberSet
        //      but at least we are not touching the simulation files :)
        Obj* obj = set.firstID();
        while (obj) 
       {
            mappe[obj->property()->name()].push_back(obj);
            obj = set.nextID(obj);
        }
        // Then we fill the dictionnary
        for (const auto &[name, group] : mappe) {
            current->objects[py::cast(name)] = group;
        }
        
    }
}

/// declare_group() : creates a python interface for an ObjGroup
template<typename Group>
auto declare_group(py::module &mod, Group group, std::string name) {
        return py::class_<Group>(mod, name.c_str())
            .def("__len__", [](const Group &v) { return v.size(); })
            .def_readwrite("prop",   &Group::prop , py::return_value_policy::reference)
            .def("__iter__", [](Group &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>())
            .def("__getitem__",[](const Group &v, size_t i) {
                     if (i >= v.size()) {
                         throw py::index_error();
                     }
                     return v[i];
                 }, py::return_value_policy::reference);
}

/// Prepares a given frame by sorting objects into object groups
Frame * prepare_frame( Simul * , int ) ;

#endif
/*
template<typename ObjProp,typename Group, typename Gmap, typename Set>
void assign_groups(Simul * sim, Gmap & mappe, std::string & name, ObjProp prop, Group group,  Set & set) {
    PropertyList plist = sim->properties.find_all(name);
    if (!plist.empty()) {
        for ( Property * i : plist )
        {
            ObjProp * fp = static_cast<ObjProp*>(i);
            mappe[fp->name()] = Group(fp);
        }
    }
    for (auto obj = set.first(); obj != set.last() ; obj = obj->next() ) {
        group[obj->property()->name()].push_back(obj);
    }
    group[set.last()->property()->name()].push_back(set.last());
}
*/
/*
/// PyObj is the python representation of an object report
struct PyObj {
    pyarray points;
    int id;
    py::dict props;
    PyObj(ObjReport* );
    PyObj() = default;
    ~PyObj() = default;
};

/// PySet is the python representation of an set report
struct PySet {
    
    py::dict props;
    py::list objects;
    PySet(SetReport* );
    PySet() = default;
    ~PySet() = default;
};

/// PySetter is the python representation of a set report
class PySetter : public py::list
{
    public:
        py::dict props;
        PySetter() = default;
        PySetter(SetReport*);
        ~PySetter() = default;
};

/// PyObjs is the python representation of a set report 
class PyObjs : public std::vector<PyObj>
{
    public:
        py::dict props;
        PyObjs() = default;
        PyObjs(SetReport*);
        ~PyObjs() = default;
};
*/

/// A class containing a vect of fibers with the same propertie
/*
class FiberGroup : public std::vector<Fiber*>
{
    public:
        FiberProp * prop;
        FiberGroup() = default;
        FiberGroup(FiberProp * p) : FiberGroup() {prop = p ;};
        ~FiberGroup() = default;
};

/// A class containing a vect of spaces with the same properties
class SpaceGroup : public std::vector<Space*>
{
    public:
        SpaceProp * prop;
        SpaceGroup() = default;
        SpaceGroup(SpaceProp * p) : SpaceGroup() {prop = p ;};
        ~SpaceGroup() = default;
};

/// A class containing a vect of spaces with the same properties
class SolidGroup : public std::vector<Solid*>
{
    public:
        SolidProp * prop;
        SolidGroup() = default;
        SolidGroup(SolidProp * p) : SolidGroup() {prop = p ;};
        ~SolidGroup() = default;
};

/// A dictionary of FiberGroups
typedef std::map<std::string, FiberGroup> Fibers;

/// A dictionary of SpaceGroups
typedef std::map<std::string, SpaceGroup> Spaces;

/// A dictionary of SolidGroups
typedef std::map<std::string, SolidGroup> Solids;
*/

/*
template<typename GMap,typename Set>
void fill_group_and_dict(Frame * current, GMap & mappe, Set & set) {
     for (auto obj = set.first(); obj != set.last() ; obj = obj->next() ) {
        mappe[obj->property()->name()].push_back(obj);
    }
    mappe[set.last()->property()->name()].push_back(set.last());
    for (const auto &[name, group] : mappe) {
        current->objects[py::cast(name)] = group;
    }
}

template<typename ObjProp,typename Group, typename Gmap>
void create_groups(Gmap & mappe, PropertyList & plist, ObjProp * prop, Group group) {
    for ( Property * i : plist )
        {
            ObjProp * fp = static_cast<ObjProp*>(i);
            mappe[fp->name()] = Group(fp);
        }
}

*/
