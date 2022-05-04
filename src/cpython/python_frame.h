#ifndef PYTHON_FRAME_H
#define PYTHON_FRAME_H

#include <fstream>
#include <sstream>
#include "sim_thread.h"
#include "stream_func.h"
//#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "organizer.h"
//#include "parser.h"
//#include "simul.h"
#include "simul_prop.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "python_utilities.h"
namespace py = pybind11;

class Simul;
class SimulProp;
class Organizer;


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
        ObjMap<Bead,BeadProp> beads;
        ObjMap<Sphere,SphereProp> spheres;
        ObjMap<Organizer,Property> organs;
        ObjMap<Space,SpaceProp> spaces;
        ObjMap<Couple,CoupleProp> couples;
        ObjMap<Single,SingleProp> singles;

        // Time of the frame
        real time;
        int index;
        int loaded;
        
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
            .def("size", &Group::size)
            .def_readwrite("prop",   &Group::prop , py::return_value_policy::reference)
            .def("__iter__", [](Group &v) {
                return py::make_iterator(v.begin(), v.end());
            }, py::keep_alive<0, 1>())
            .def("__getitem__",[](const Group &v, size_t i) {
                int s = v.size();
                if (i<0) {i+=s;} // Python time negative indexing
				if (i >= s or i<0) {
                         throw py::index_error();
                     }
                     return v[i];
                 }, py::return_value_policy::reference);
}

/// Prepares a given frame by sorting objects into object groups
Frame * make_frame( Simul * sim) 
{   
    std::vector<std::string> categories = std::vector<std::string>{"aster","nucleus","bundle","fake"};
    //extern std::vector<std::string>  categories;
    Frame * current = new Frame;
    current->simul = sim;
    
    distribute_objects(sim,current, current->fibers, sim->fibers, std::string("fiber") ) ;
    distribute_objects(sim,current, current->solids, sim->solids, std::string("solid") ) ;
    distribute_objects(sim,current, current->spaces, sim->spaces, std::string("space") ) ;
    distribute_objects(sim,current, current->beads, sim->beads, std::string("bead") ) ;
    distribute_objects(sim,current, current->spheres, sim->spheres, std::string("sphere") ) ;
    // For organizer, the we have to check the different categories
    for (auto categ : categories) {
        distribute_objects(sim,current, current->organs, sim->organizers, std::string(categ) ) ;
    }
    // for couple and single we need to use firstID, nextID
    distribute_objects_wID(sim,current, current->couples, sim->couples, std::string("couple") ) ;
    distribute_objects_wID(sim,current, current->singles, sim->singles, std::string("single") ) ;
    
    current->time = sim->time();
    //current->index = frame;
    current->loaded = 1;

    return current;
    
}

/// Prepares a given frame and attibutes an index
Frame * make_frame_index( Simul * sim, int i)  {
	Frame * current = make_frame(sim);
	current->index = i;
	return current;
}


#endif
