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
    
    // This should probably infor fiber_modules.h ?
    py::class_<FiberProp,Property>(m, "FiberProp")
        .def_readwrite("segmentation", &FiberProp::segmentation)
        .def_readwrite("rigidity", &FiberProp::rigidity)
        .def_readwrite("min_length", &FiberProp::min_length)
        .def_readwrite("max_length", &FiberProp::max_length)
        .def_readwrite("total_polymer", &FiberProp::total_polymer)
        .def_readwrite("persistent", &FiberProp::persistent)
        .def_readwrite("viscosity", &FiberProp::viscosity)
        .def_readwrite("drag_radius", &FiberProp::drag_radius)
        .def_readwrite("drag_length", &FiberProp::drag_length)
        .def_readwrite("drag_model", &FiberProp::drag_model)
        .def_readwrite("drag_gap", &FiberProp::drag_gap)
        .def_readwrite("binding_key", &FiberProp::binding_key)
        .def_readwrite("lattice", &FiberProp::lattice)
        .def_readwrite("lattice_unit", &FiberProp::lattice_unit)
        .def_readwrite("confine", &FiberProp::confine)
        .def_readwrite("confine_space", &FiberProp::confine_space)
        .def_readwrite("confine_stiffness", &FiberProp::confine_stiffness)
        .def_readwrite("steric", &FiberProp::steric)
        .def_readwrite("steric_radius", &FiberProp::steric_radius)
        .def_readwrite("steric_range", &FiberProp::steric_range)
        .def_readwrite("glue", &FiberProp::glue)
        .def_readwrite("glue_single", &FiberProp::glue_single)
        .def_readwrite("activity", &FiberProp::activity)
        .def_readwrite("display_fresh", &FiberProp::display_fresh)
        .def_readwrite("display", &FiberProp::display);
        
    py::class_<HandProp,Property>(m, "HandProp");
    py::class_<CoupleProp,Property>(m, "CoupleProp");
    py::class_<SingleProp,Property>(m, "SingleProp");
    py::class_<SpaceProp,Property>(m, "SpaceProp");
}

