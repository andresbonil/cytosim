#include "space.h"
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
}

