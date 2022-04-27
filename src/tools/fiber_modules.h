#include "fiber.h"
#include "classic_fiber.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Fiber;
class Object;

/// a utility to enrich the cytosim python module
void load_fiber_classes(py::module_ &m) {
    /// Python interface to Fiber
    // @TODO : add the methods from mecafil, chain, mecable... ?
    py::class_<Fiber,Object>(m, "Fiber")
        .def("points",  [](Fiber * fib) {return get_obj_points(fib);})
        .def("nbPoints",  [](Fiber * fib) {return fib->nbPoints();})
        .def("cutM",  [](Fiber * fib, real len) {return fib->cutM(len);})
        .def("cutP",  [](Fiber * fib, real len) {return fib->cutP(len);})
        .def("sever",  [](Fiber * fib, real a, int p, int m) {return fib->sever(a, p, m);})
        .def("severNow",  [](Fiber * fib, real a) {return fib->severNow(a);})
        .def("severP",  [](Fiber * fib, real a) {return fib->severP(a);})
        .def("nowSever",  [](Fiber * fib) {return fib->severNow();})
        .def("join",  [](Fiber * fib, Fiber * fob) {return fib->join(fob);})
        .def("stateM",  [](const Fiber * fib) {return py::cast(fib->dynamicStateM());})
        .def("stateP",  [](const Fiber * fib) {return py::cast(fib->dynamicStateP());})
        .def("setStateM",  [](Fiber * fib, int stat) {return fib->setDynamicStateM(stat);})
        .def("setStateP",  [](Fiber * fib, int stat) {return fib->setDynamicStateP(stat);})
        .def("setDynamicState",  [](Fiber * fib, int end, int stat) {return fib->setDynamicState(static_cast<FiberEnd>(end),stat);})
        .def("freshAssembly",  [](Fiber * fib, int end) {return fib->freshAssembly(static_cast<FiberEnd>(end));})
        .def("nbHands",  [](Fiber * fib) {return fib->nbHands();})
        .def("nbHandsInRange",  [](Fiber * fib, real amin, real amax, int end) {return fib->nbHandsInRange(amin, amax, static_cast<FiberEnd>(end));})
        .def("__next__", [](const Fiber * fib) {return fib->next();});
        
        
    py::class_<ClassicFiber,Fiber>(m, "ClassicFiber")
        .def("freshAssemblyM",  [](ClassicFiber * fib) {return fib->freshAssemblyM();})
        .def("freshAssemblyP",  [](ClassicFiber * fib) {return fib->freshAssemblyP();});
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
         
    /// Python interface to FiberProp
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
        
}

