#include "solid.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Solid;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_solid_classes(py::module_ &m) {
     /// Python interface to Fiber
    py::class_<Solid,Object>(m, "Fiber")
        .def("position", [](const Solid * sol) {return to_numpy(sol->position());})
        .def("id",  [](const Solid * sol) {return sol->identity();})
        .def("info",  [](const Solid * sol) {return to_dict(sol->info());});
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

