#include "solid.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Solid;
class Object;

/// a utility to enrich the cytosim python module
void load_solid_classes(py::module_ &m) {
     /// Python interface to Solid
    py::class_<Solid,Object>(m, "Solid")
        .def("position", [](const Solid * sol) {return to_numpy(sol->position());})
        .def("points",  [](Solid * sol) {return get_obj_points(sol);})
        .def("dragCoefficient",  [](Solid * sol) {return sol->dragCoefficient() ;})
        .def("addSphere",  [](Solid * sol, pyarray pts, real radius) {return sol->addSphere(to_vector(pts), radius) ;})
        .def("setRadius",  [](Solid * sol, int i, real radius) {return sol->radius((unsigned) i, radius) ;})
        .def("radius",  [](Solid * sol, int i) {return sol->radius((unsigned) i) ;})
        .def("centroid", [](const Solid * sol) {return to_numpy(sol->centroid());})
        .def("next",  [](Solid * sol) {return sol->next() ;})
        .def("prev",  [](Solid * sol) {return sol->prev() ;})
        .def("tag",  [](Solid * sol) {return sol->tag() ;})
        .def("property",  [](Solid * sol) {return sol->property() ;})
        .def("nbPoints", [](const Solid * sol) {return sol->nbPoints();});
        /**
         * 
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

