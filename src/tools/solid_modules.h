#include "solid.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
#include "objecter_python.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Solid;
class Object;
//class Property;

real_array * get_points(Solid* sol) {
    int_vect sizes = {(int)sol->nbPoints(), (int)DIM};
    int_vect strides = {DIM*sizeof(real), sizeof(real)};
    return new real_array{sol->data(), sizes, strides};;
}

/// a utility to enrich the cytosim python module
void load_solid_classes(py::module_ &m) {
     /// Python interface to Solid
    py::class_<Solid,Object>(m, "Solid")
        .def("position", [](const Solid * sol) {return to_numpy(sol->position());})
        .def("points",  [](Solid * sol) {return to_numpy(get_points(sol));})
        .def("nbPoints", [](const Solid * sol) {return sol->nbPoints();});
        
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
}

