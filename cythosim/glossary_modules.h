#include "aster.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Aster;
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

    py::class_<Aster,Organizer>(m, "Aster")
        .def("solid",  [](const Aster * org) {return org->solid() ;})
        .def("position", [](const Aster * org) {return to_numpy(org->position());})
        .def("nbFibers",  [](const Aster * org) {return org->nbFibers() ;})
        .def("fiber",  [](const Aster * org, int n) {return org->fiber((size_t)n) ;})
        .def("posLink1",  [](const Aster * org, int n) {return to_numpy(org->posLink1((size_t)n)) ;})
        .def("posLink2",  [](const Aster * org, int n) {return to_numpy(org->posLink2((size_t)n)) ;})
        .def("posFiber1",  [](const Aster * org, int n) {return to_numpy(org->posFiber1((size_t)n)) ;})
        .def("posFiber2",  [](const Aster * org, int n) {return to_numpy(org->posFiber2((size_t)n)) ;})
        .def("getLink1",  [](const Aster * org, int n)
            {Vector V,W; 
            org->getLink1((size_t)n,V,W);
                return std::vector<pyarray>{to_numpy(V),to_numpy(W)};
                })
        .def("getLink2",  [](const Aster * org, int n)
            {Vector V,W; 
            org->getLink2((size_t)n,V,W);
                return std::vector<pyarray>{to_numpy(V),to_numpy(W)};
                })
        .def("getLink",  [](const Aster * org, int n)
            {Vector V,W; 
            org->getLink((size_t)n,V,W);
                return std::vector<pyarray>{to_numpy(V),to_numpy(W)};
                })
        .def("property",  [](const Aster * org) {return org->property() ;}); 
}

