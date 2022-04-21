#include "hand.h"
#include "fiber.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Hand;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_hand_classes(py::module_ &m) {
     /// Python interface to couple
    py::class_<Hand>(m, "Hand")
        .def("fiber",  [](Hand * h) {return h->fiber();})
        .def("property",  [](Hand * h) {return h->property();})
        .def("relocate",  [](Hand * h, Fiber * fib) {return h->relocate(fib);})
        .def("relocate_to",  [](Hand * h, Fiber * fib, real a) {return h->relocate(fib,a);})
        .def("moveToEnd",  [](Hand * h, int end) {return h->moveToEnd(static_cast<FiberEnd>(end));})
        .def("locate",  [](Hand * h, Fiber * fib, real a) {return h->locate(fib,a);})
        .def("attach",  [](Hand * h, Fiber * fib, real a,  int end) {return h->attach(fib, a, static_cast<FiberEnd>(end));});
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

