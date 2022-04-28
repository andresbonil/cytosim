#include "fiber_prop.h"
#include "simul_prop.h"
#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Property;

/**
 * @brief 
 * @TODO : manage to have objectSet ! Now missing operator = ?////

 */

/// a utility to enrich the cytosim python module
auto load_simul_classes(py::module_ &m) {
    /// Python interface to default property
    py::class_<Property>(m, "Prop")
        .def("name", &Property::name)
        .def("complete",  [](Property * prop, Simul * sim) {return prop->complete(*sim);});
    
    py::class_<ObjectSet>(m, "ObjectSet");
    py::class_<SpaceSet,ObjectSet>(m, "SpaceSet");
    
    auto pysim = py::class_<Simul>(m, "Simul")
        .def_readwrite("prop",   &Simul::prop , py::return_value_policy::reference)
        .def_readwrite("properties",   &Simul::properties , py::return_value_policy::reference)
        //.def_readwrite("fields",   &Simul::fields , py::return_value_policy::reference)
        //.def_readwrite("fibers",   &Simul::fibers , py::return_value_policy::reference)
        //.def_readwrite("spheres",   &Simul::spheres , py::return_value_policy::reference)
        //.def_readwrite("beads",   &Simul::beads , py::return_value_policy::reference)
        //.def_readwrite("solids",   &Simul::solids , py::return_value_policy::reference)
        //.def_readwrite("couples",   &Simul::couples , py::return_value_policy::reference)
        //.def_readwrite("singles",   &Simul::singles , py::return_value_policy::reference)
        //.def_readwrite("organizers",   &Simul::organizers , py::return_value_policy::reference)
        .def("remove",  [](Simul * sim, Object* obj) {return sim->remove(obj);})
        .def("erase",  [](Simul * sim, Object* obj) {return sim->erase(obj);})
        .def("nuke",  [](Simul * sim) {return sim->erase();})
        .def("time",  [](Simul * sim) {return sim->time();})
        .def("time_step",  [](Simul * sim) {return sim->time_step();})
        .def("step",  [](Simul * sim) {return sim->step();})
        .def("solve",  [](Simul * sim) {return sim->solve();})
        .def("solve_auto",  [](Simul * sim) {return sim->solve_auto();})
        //.def("dump",  [](Simul * sim, std::string s) {return sim->dump( &s[0]);})
        //.def("saveSystem",  [](Simul * sim, char s) {return sim->saveSystem((char) s);})
        .def("evaluate",  [](Simul * sim, std::string s) {return sim->evaluate(s);} , py::return_value_policy::reference)
        .def("toMecable",  [](Simul * sim, Object* o) {return sim->toMecable(o);} , py::return_value_policy::reference)
        .def("findMecable",  [](Simul * sim, std::string s) {return sim->findMecable(s);} , py::return_value_policy::reference)
        .def("findSpace",  [](Simul * sim, std::string s) {return sim->findSpace(s);} , py::return_value_policy::reference)
        .def("rename",  [](Simul * sim, std::string s) {return sim->rename(s);} , py::return_value_policy::reference)
        .def("isCategory",  [](Simul * sim, std::string s) {return sim->isCategory(s);} , py::return_value_policy::reference)
        .def("findProperty",  [](Simul * sim, std::string s) {return sim->findProperty(s);} , py::return_value_policy::reference);
        
        
        
    
    /// Python interface to simulProp
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

    return pysim;
}

