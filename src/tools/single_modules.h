#include "single.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Single;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_single_classes(py::module_ &m) {
     /// Python interface to single
    py::class_<Single,Object>(m, "Single")
        .def("state",  [](Single * s) {return s->state();});
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

