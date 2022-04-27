#include "objecter.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Object;

/// a utility to enrich the cytosim python module
void load_object_classes(py::module_ &m) {
     /// Python interface to Organizer
    py::class_<Object>(m, "Object")
        .def("reference",  [](Object * obj) {return obj->reference() ;})
        .def("property",  [](Object * obj) {return obj->property() ;})
        .def("position", [](const Object * obj) {return to_numpy(obj->position());})
        .def("next",  [](Object * obj) {return obj->next() ;})
        .def("prev",  [](Object * obj) {return obj->prev() ;})
        .def("id",  [](const Object * obj) {return obj->identity();})
        .def("points", [](const Object * obj) {return pyarray();});
        
}
