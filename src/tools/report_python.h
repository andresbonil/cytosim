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
#include "python_utilities.h"

namespace py = pybind11;


class Simul;
class SimulProp;

/*
class ObjWrapper {
    public:
        Object * obj;
        ObjWrapper(Object * o) { obj = o; };
        ~ObjWrapper() {};
};
*/

/// A class containing a vect of fibers with the same properties
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

/// A time frame ; basicaly 
class Frame 
{
    //py::dict objects;
    public:
        Fibers fibers;
        Spaces spaces;
        Solids solids;
        int time;
        py::dict objects;
        Frame() = default;
        ~Frame() = default;
};

//Frame & prepare_frame(int ) ;
//Frame & prepare_frame( Simul * , int ) ;
Frame * prepare_frame( Simul * , int ) ;

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