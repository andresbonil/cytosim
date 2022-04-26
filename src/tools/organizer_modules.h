#include "organizer.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Organizer;
class Object;

/// a utility to enrich the cytosim python module
void load_organizer_classes(py::module_ &m) {
     /// Python interface to Organizer
    py::class_<Organizer,Object>(m, "Organizer")
        .def("nbOrganized",  [](const Organizer * org) {return org->nbOrganized() ;})
        .def("position", [](const Organizer * org) {return to_numpy(org->position());})
        .def("positionP", [](const Organizer * org, unsigned i) {return to_numpy(org->positionP(i));})
        .def("dragCoefficient",  [](Organizer * org) {return org->dragCoefficient() ;})
        .def("next",  [](Organizer * org) {return org->next() ;})
        .def("prev",  [](Organizer * org) {return org->prev() ;});


}

