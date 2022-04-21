#include "couple.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Couple;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_couple_classes(py::module_ &m) {
     /// Python interface to couple
    py::class_<Couple,Object>(m, "Couple")
        .def("state",  [](Couple * s) {return s->state();});
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

