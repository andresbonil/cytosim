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
    py::class_<Mecable,Object>(m, "Mecable")
        .def("nbPoints", [](const Mecable * mec) {return mec->nbPoints();})
        .def("allocated", [](const Mecable * mec) {return mec->allocated();})
        .def("posPoint", [](const Mecable * mec, unsigned p) {return to_numpy(mec->posPoint(p));})
        .def("points",  [](const Mecable * mec) {return get_obj_points(mec);})
        .def("setPoint", []( Mecable * mec, unsigned p, pyarray x) {return mec->setPoint(p,to_vector(x));})
        .def("movePoint", []( Mecable * mec, unsigned p, pyarray x) {return mec->movePoint(p,to_vector(x));})
        .def("addPoint", []( Mecable * mec, pyarray x) {return mec->addPoint(to_vector(x));})
        .def("removePoints", []( Mecable * mec, unsigned p, unsigned q) {return mec->removePoints(p,q);})
        .def("clearPoints", []( Mecable * mec) {return mec->clearPoints();})
        .def("shiftPoints", []( Mecable * mec, unsigned p, unsigned q) {return mec->shiftPoints(p,q);})
        .def("truncateM", []( Mecable * mec, unsigned p) {return mec->truncateM(p);})
        .def("truncateP", []( Mecable * mec, unsigned p) {return mec->truncateP(p);})
        .def("calculateMomentum1", [](Mecable * mec, bool sub) 
            {Vector V,W; 
            mec->calculateMomentum(V,W,sub);
            return std::vector<pyarray>{to_numpy(V),to_numpy(W)};
            })
        .def("dragCoefficient", [](const Mecable * mec) {return mec->dragCoefficient();})
        .def("netForce", [](const Mecable * mec, unsigned p) {return to_numpy(mec->netForce(p));})
        .def("position", [](const Mecable * mec) {return to_numpy(mec->position());})
        .def("mobile", [](const Mecable * mec) {return mec->mobile();})
        .def("translate", []( Mecable * mec, pyarray x) {return mec->translate(to_vector(x));})
        .def("allInside", []( Mecable * mec, Space * x) {return mec->allInside(x);});
        
    py::class_<Sphere,Mecable>(m, "Sphere");
    
    py::class_<Bead,Mecable>(m, "Bead");
    
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
            @TODO : ADD SPECIALIZED SOLID CLASSES
         */
}

