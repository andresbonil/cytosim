#include "hand.h"
#include "fiber.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Hand;
class FiberSite;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_hand_classes(py::module_ &m) {
     /// Python interface to Fibersite
     /*
      Now unused because of a memory bug
    */
    py::class_<FiberSite>(m,"FiberSite")
        .def("moveTo",  [](FiberSite * h, real a) {return h->moveTo(a);})
        .def("relocateM",  [](FiberSite * h) {return h->relocateM();})
        .def("relocateP",  [](FiberSite * h) {return h->relocateP();})
        .def("unattached",  [](FiberSite * h) {return h->unattached();})
        .def("attached",  [](FiberSite * h) {return h->attached();})
        .def("update",  &FiberSite::update)
        .def("interpolation",  &FiberSite::interpolation, py::return_value_policy::reference)
        .def("fiber",  [](FiberSite * h) {return h->fiber();}, py::return_value_policy::reference)
        .def("position",  [](FiberSite * h) {return to_numpy(h->pos());})
        .def("posHand",  [](FiberSite * h) {return to_numpy(h->posHand());})
        .def("direction",  [](FiberSite * h) {return to_numpy(h->dir());}) // direction because dir has a python meaning
        .def("dirFiber",  [](FiberSite * h) {return to_numpy(h->dirFiber());})
        .def("abscissa",  [](FiberSite * h) {return h->abscissa();})
        .def("abscissaFromM",  [](FiberSite * h) {return h->abscissaFromM();})
        .def("abscissaFromP",  [](FiberSite * h) {return h->abscissaFromP();})
        .def("abscissaFrom",  [](FiberSite * h, int end) {return h->abscissaFrom(static_cast<FiberEnd>(end));})
        .def("nearestEnd",  [](FiberSite * h) {return static_cast<int>(h->nearestEnd());})
        .def("distanceToEnd",  [](FiberSite * h, int end) {return h->distanceToEnd(static_cast<FiberEnd>(end));});
    
        
    py::class_<Hand,FiberSite>(m, "Hand")
        .def("property",  [](Hand * h) {return h->property();})
        .def("prop",  [](Hand * h) {return h->prop;})
        .def("relocate",  [](Hand * h, Fiber * fib) {return h->relocate(fib);})
        .def("relocateTo",  [](Hand * h, Fiber * fib, real a) {return h->relocate(fib,a);})
        .def("moveToEnd",  [](Hand * h, int end) {return h->moveToEnd(static_cast<FiberEnd>(end));})
        .def("locate",  [](Hand * h, Fiber * fib, real a) {return h->locate(fib,a);})
        .def("attach",  [](Hand * h, Fiber * fib, real a,  int end) {return h->attach(fib, a, static_cast<FiberEnd>(end));})
        .def("detach",  [](Hand * h) {return h->detach();})
        .def("attachEnd",  [](Hand * h, Fiber * fib, int end) {return h->attachEnd(fib, static_cast<FiberEnd>(end));})
        .def("attachTo",  [](Hand * h, Fiber * fib, real a,  int end) {return h->attachTo(fib, a, static_cast<FiberEnd>(end));})
        .def("otherHand",  [](Hand * h) {return h->otherHand();})
        .def("otherPosition",  [](Hand * h) {return to_numpy(h->otherPosition());})
        .def("linkStiffness",  [](Hand * h) {return h->linkStiffness();})
        /* 
         These should be in FiberSite, but there is a memory bug 
        */
        .def("moveTo",  [](Hand * h, real a) {return h->moveTo(a);})
        .def("relocateM",  [](Hand * h) {return h->relocateM();})
        .def("relocateP",  [](Hand * h) {return h->relocateP();})
        .def("unattached",  [](Hand * h) {return h->unattached();})
        .def("attached",  [](Hand * h) {return h->attached();})
        .def("update",  &Hand::update)
        .def("interpolation",  &Hand::interpolation, py::return_value_policy::reference)
        .def("fiber",  [](Hand * h) {return h->fiber();}, py::return_value_policy::reference)
        .def("position",  [](Hand * h) {return to_numpy(h->pos());})
        .def("posHand",  [](Hand * h) {return to_numpy(h->posHand());})
        .def("direction",  [](Hand * h) {return to_numpy(h->dir());}) // direction because dir has a python meaning
        .def("dirFiber",  [](Hand * h) {return to_numpy(h->dirFiber());})
        .def("abscissa",  [](Hand * h) {return h->abscissa();})
        .def("abscissaFromM",  [](Hand * h) {return h->abscissaFromM();})
        .def("abscissaFromP",  [](Hand * h) {return h->abscissaFromP();})
        .def("abscissaFrom",  [](Hand * h, int end) {return h->abscissaFrom(static_cast<FiberEnd>(end));})
        .def("nearestEnd",  [](Hand * h) {return static_cast<int>(h->nearestEnd());})
        .def("distanceToEnd",  [](Hand * h, int end) {return h->distanceToEnd(static_cast<FiberEnd>(end));});
         
        
         /**
            @TODO : ADD SPECIALIZED HAND CLASSES
         */
    py::class_<HandProp,Property>(m, "HandProp")
        .def_readwrite("binding_rate", &HandProp::binding_rate)
        .def_readwrite("binding_range", &HandProp::binding_range)
        .def_readwrite("binding_key", &HandProp::binding_key)
        .def_readwrite("unbinding_rate", &HandProp::unbinding_rate)
        .def_readwrite("unbinding_force", &HandProp::unbinding_force)
        .def_readwrite("bind_also_end", &HandProp::bind_also_end)
        .def_readwrite("bind_end_range", &HandProp::bind_end_range)
        .def_readwrite("hold_growing_end", &HandProp::hold_growing_end)
        .def_readwrite("hold_shrinking_end", &HandProp::hold_shrinking_end)
        .def_readwrite("activity", &HandProp::activity)
        .def_readwrite("display", &HandProp::display);
        
}

