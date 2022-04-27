#include "fiber_prop.h"
#include "simul_prop.h"
#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class FiberProp;
class Property;

/// a utility to enrich the cytosim python module
void load_prop_classes(py::module_ &m) {
    py::class_<Property>(m, "Prop")
        .def("name", &Property::name); // prop.name() outputs the name
    
    py::class_<SimulProp,Property>(m, "SimulProp")
        .def_readwrite("time", &SimulProp::time)
        .def_readwrite("time_step", &SimulProp::time_step)
        .def_readwrite("viscosity", &SimulProp::viscosity)
        .def_readwrite("kT", &SimulProp::kT)
        .def_readwrite("tolerance", &SimulProp::tolerance)
        .def_readwrite("acceptable_prop", &SimulProp::acceptable_prob)
        .def_readwrite("steric", &SimulProp::steric)
        //.def_readwrite("steric_stiffness_push", &SimulProp::steric_stiffness_push) <- issue with real[]
        //.def_readwrite("steric_stiffness_pull", &SimulProp::steric_stiffness_pull)
        .def_readwrite("steric_max_range", &SimulProp::steric_max_range)
        .def_readwrite("verbose", &SimulProp::verbose)
        .def_readwrite("config_file", &SimulProp::config_file)
        .def_readwrite("trajectory_file", &SimulProp::trajectory_file)
        .def_readwrite("clear_trajectory", &SimulProp::clear_trajectory)
        .def_readwrite("skip_free_couple", &SimulProp::skip_free_couple)
        .def_readwrite("display_fresh", &SimulProp::display_fresh)
        .def_readwrite("display", &SimulProp::display);
    

}

