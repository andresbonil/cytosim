#include "fiber.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Fiber;
class Object;
//class Property;

/// a utility to enrich the cytosim python module
void load_fiber_classes(py::module_ &m) {
     /// Python interface to Fiber
    py::class_<Fiber,Object>(m, "Fiber")
        .def("stateM",  [](const Fiber * fib) {return py::cast(fib->dynamicStateM());})
        .def("stateP",  [](const Fiber * fib) {return py::cast(fib->dynamicStateP());})
        .def("setStateM",  [](Fiber * fib, int stat) {return fib->setDynamicStateM(stat);})
        .def("setStateP",  [](Fiber * fib, int stat) {return fib->setDynamicStateP(stat);})
        .def("__next__", [](const Fiber * fib) {return fib->next();});
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

