#include "interpolation.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Interpolation;

/// a utility to enrich the cytosim python module
void load_point_classes(py::module_ &m) {
    py::class_<Mecapoint>(m, "Mecapoint")
        .def("set", &Mecapoint::set)
        .def("mecable", &Mecapoint::mecable, py::return_value_policy::reference)
        .def("valid", &Mecapoint::valid)
        .def("position", [](const Mecapoint * pol) {return to_numpy(pol->pos());})
        .def("pos", [](const Mecapoint * pol) {return to_numpy(pol->pos());})
        .def("overlapping", &Mecapoint::overlapping)
        .def("near", &Mecapoint::near);
        
    py::class_<Interpolation>(m, "Interpolation")
        .def("clear", &Interpolation::clear)
        .def("set", &Interpolation::set)
        .def("valid", &Interpolation::valid)
        .def("mecable", &Interpolation::mecable, py::return_value_policy::reference)
        .def("exact1", &Interpolation::exact1, py::return_value_policy::reference)
        .def("exact2", &Interpolation::exact2, py::return_value_policy::reference)
        .def("coef1", &Interpolation::coef1)
        .def("coef2", &Interpolation::coef2)
        .def("coef", &Interpolation::coef)
        .def("position", [](const Interpolation * pol) {return to_numpy(pol->pos());})
        .def("pos", [](const Interpolation * pol) {return to_numpy(pol->pos());})
        .def("pos1", [](const Interpolation * pol) {return to_numpy(pol->pos1());})
        .def("pos2", [](const Interpolation * pol) {return to_numpy(pol->pos2());})
        .def("diff", [](const Interpolation * pol) {return to_numpy(pol->diff());})
        .def("len", &Interpolation::len)
        .def("lenSqr", &Interpolation::lenSqr)
        .def("dir", [](const Interpolation * pol) {return to_numpy(pol->dir());})
        .def("inside", &Interpolation::inside)
        .def("overlapping", [](const Interpolation * pol, const Interpolation & pal) {return pol->overlapping(pal);})
        .def("overlapPoint", [](const Interpolation * pol, const Mecapoint & pal) {return pol->overlapping(pal);});

}
