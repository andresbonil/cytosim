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
     /// Python interface to couple
    py::class_<FiberSite>(m,"FiberSite")
        .def("moveTo",  [](FiberSite * h, real a) {return h->moveTo(a);})
        .def("relocateM",  [](FiberSite * h) {return h->relocateM();})
        .def("relocateP",  [](FiberSite * h) {return h->relocateP();})
        .def("unattached",  [](FiberSite * h) {return h->unattached();})
        .def("attached",  [](FiberSite * h) {return h->attached();})
        .def("fiber",  [](FiberSite * h) {return h->fiber();})
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
        .def("linkStiffness",  [](Hand * h) {return h->linkStiffness();});
        
        //.def("attachEnd2",  [](Couple * s, Fiber * fib, int end) {return s->attachEnd2(fib, static_cast<FiberEnd>(end));})
        //.def("moveToEnd1",  [](Couple * s,int end) {return s->moveToEnd1(static_cast<FiberEnd>(end));})
        //.def("moveToEnd2",  [](Couple * s,int end) {return s->moveToEnd2(static_cast<FiberEnd>(end));})
        //.def("fiber1",  [](Couple * s) {return s->fiber1();})
        //.def("fiber2",  [](Couple * s) {return s->fiber2();})
        //.def("abcissa",  [](Couple * s) {return to_numpy(s->posFree());})
        
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

