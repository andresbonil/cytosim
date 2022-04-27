#include "space.h"
#include "common.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Space;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_space_classes(py::module_ &m) {
     /// Python interface to space
    py::class_<Space,Object>(m, "Space");
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
    py::class_<SpaceProp,Property>(m, "SpaceProp");
            
    py::enum_<Confinement>(m,"Confinement")
        .value("CONFINE_OFF", CONFINE_OFF)
        .value("CONFINE_INSIDE", CONFINE_INSIDE)
        .value("CONFINE_OUTSIDE", CONFINE_OUTSIDE)
        .value("CONFINE_ON", CONFINE_ON)
        .value("CONFINE_ALL_INSIDE", CONFINE_ALL_INSIDE)
        .value("CONFINE_ALL_OUTSIDE", CONFINE_ALL_OUTSIDE)
        .value("CONFINE_PLUS_END", CONFINE_PLUS_END)
        .value("CONFINE_MINUS_END", CONFINE_MINUS_END)
        .value("CONFINE_BOTH_ENDS", CONFINE_BOTH_ENDS)
        .value("CONFINE_PLUS_OUT", CONFINE_PLUS_OUT)
        .value("CONFINE_POINT", CONFINE_POINT)
        .value("CONFINE_RANGE", CONFINE_RANGE)
        .export_values();

}

