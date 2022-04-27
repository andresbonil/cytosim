#include "couple.h"
#include "fiber.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Couple;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_couple_classes(py::module_ &m) {
     /// Python interface to couple
    py::class_<Couple,Object>(m, "Couple")
        .def("position",  [](Couple * s) {return to_numpy(s->position());})
        .def("active",  [](Couple * s) {return s->active();})
        .def("stiffness",  [](Couple * s) {return s->stiffness();})
        .def("force",  [](Couple * s) {return to_numpy(s->force());})
        .def("cosAngle",  [](Couple * s) {return s->cosAngle();})
        .def("sidePos",  [](Couple * s) {return to_numpy(s->sidePos());})
        .def("posFree",  [](Couple * s) {return to_numpy(s->posFree());})
        // all this definitely should be in the interface to hand
        //.def("attachEnd1",  [](Couple * s, Fiber * fib, int end) {return s->attachEnd1(fib, static_cast<FiberEnd>(end));})
        //.def("attachEnd2",  [](Couple * s, Fiber * fib, int end) {return s->attachEnd2(fib, static_cast<FiberEnd>(end));})
        //.def("moveToEnd1",  [](Couple * s,int end) {return s->moveToEnd1(static_cast<FiberEnd>(end));})
        //.def("moveToEnd2",  [](Couple * s,int end) {return s->moveToEnd2(static_cast<FiberEnd>(end));})
        //.def("fiber1",  [](Couple * s) {return s->fiber1();})
        //.def("fiber2",  [](Couple * s) {return s->fiber2();})
        //.def("abcissa",  [](Couple * s) {return to_numpy(s->posFree());})
        .def("hand1",  [](Couple * s) {return s->hand1();}, py::return_value_policy::reference)
        .def("hand2",  [](Couple * s) {return s->hand2();}, py::return_value_policy::reference)
        .def("hand",  [](Couple * s, int i) {
            if (i==0) {return s->hand1();} else {return s->hand2();} ;}
            , py::return_value_policy::reference)
        .def("state",  [](Couple * s) {return s->state();})
        .def("__getitem__",[](const Couple *s, int i) { // We can call couple[0]  to get the first hand ! thus couple[0].attachEnd(...) is available
            if (i==0) {return s->hand1();} else {return s->hand2();} ;}
            , py::return_value_policy::reference);
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
         
        py::class_<CoupleProp,Property>(m, "CoupleProp")
            .def_readwrite("hand1", &CoupleProp::hand1)
            .def_readwrite("hand2", &CoupleProp::hand2)
            .def_readwrite("stiffness", &CoupleProp::stiffness)
            .def_readwrite("length", &CoupleProp::length)
            .def_readwrite("diffusion", &CoupleProp::diffusion)
            .def_readwrite("fast_diffusion", &CoupleProp::fast_diffusion)
            .def_readwrite("stiff", &CoupleProp::stiff)
            .def_readwrite("specificity", &CoupleProp::specificity)
            .def_readwrite("confine", &CoupleProp::confine)
            .def_readwrite("confine_space", &CoupleProp::confine_space)
            .def_readwrite("activity", &CoupleProp::activity)
            .def_readwrite("hand1_prop", &CoupleProp::hand1_prop)
            .def_readwrite("hand2_prop", &CoupleProp::hand2_prop);
            
            
}

