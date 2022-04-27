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
    py::class_<SingleProp,Property>(m, "SingleProp")
        .def_readwrite("hand", &SingleProp::hand)
        .def_readwrite("stiffness", &SingleProp::stiffness)
        .def_readwrite("length", &SingleProp::length)
        .def_readwrite("diffusion", &SingleProp::diffusion)
        .def_readwrite("fast_diffusion", &SingleProp::fast_diffusion)
        .def_readwrite("confine", &SingleProp::confine)
        .def_readwrite("confine_space", &SingleProp::confine_space)
        .def_readwrite("activity", &SingleProp::activity)
        .def_readwrite("hand_prop", &SingleProp::hand_prop);
}

