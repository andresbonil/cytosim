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
        .def("data",   [](Mecable * mec) {return get_pointsarray(mec);})
        .def("nbPoints", &Mecable::nbPoints)
        .def("allocated", &Mecable::allocated)
        .def("points",  [](Mecable * mec) {return get_obj_points(mec);})
        .def("posPoint",  [](Mecable * mec,int p) {return to_numpy(mec->posPoint(p));})
        .def("setPoint",  [](Mecable * mec, int i, pyarray vec) 
            {   Vector p = to_vector(vec);
                return mec->setPoint(i,p);})
        .def("setPoint",  [](Mecable * mec, int i, pyarray vec) 
            {   Vector p = to_vector(vec);
                return mec->movePoint(i,p);})
        .def("addPoint",  [](Mecable * mec, pyarray vec) 
            {   Vector p = to_vector(vec);
                return mec->addPoint(p);})
        .def("removePoints", &Mecable::removePoints)
        .def("clearPoints", &Mecable::clearPoints)
        .def("shiftPoints", &Mecable::shiftPoints)
        .def("truncateM", &Mecable::truncateM)
        .def("truncateP", &Mecable::truncateP)
        .def("calculateMomentum",  [](Mecable * mec, bool sub)
            {Vector V,W; 
            mec->calculateMomentum(V,W,sub);
                return std::vector<pyarray>{to_numpy(V),to_numpy(W)};
                })
        .def("netForce",  [](Mecable * mec,int p) {return to_numpy(mec->netForce(p));})
        .def("position",  [](Mecable * mec) {return to_numpy(mec->position());})
        .def("translate",  [](Mecable * mec, pyarray vec) 
            {   Vector p = to_vector(vec);
                return mec->translate(p);})
        .def("allInside", &Mecable::allInside);
        
    py::class_<Sphere,Mecable>(m, "Sphere")
        .def("position",  [](Sphere * bed) {return to_numpy(bed->position());})
        .def("pos",  [](Sphere * bed) {return to_numpy(bed->position());})
        .def("radius", &Sphere::radius)
        .def("resize", &Sphere::resize)
        .def("reshape", &Sphere::reshape)
        .def("orthogonalize", &Sphere::orthogonalize)
        .def("addSurfacePoint", [](Sphere & bed, pyarray pos) {bed.addSurfacePoint(to_vector(pos));})
        .def("nbbSurfacePoints", &Sphere::nbSurfacePoints)
        .def("dragCoefficient", &Sphere::dragCoefficient)
        .def("next", &Sphere::next, py::return_value_policy::reference)
        .def("prev", &Sphere::prev, py::return_value_policy::reference)
        .def("toSphere",  [](Object * obj) {return Sphere::toSphere(obj);},  py::return_value_policy::reference);
        
    
    py::class_<Bead,Mecable>(m, "Bead")
        .def("position",  [](Bead * bed) {return to_numpy(bed->position());})
        .def("pos",  [](Bead * bed) {return to_numpy(bed->position());})
        .def("setPosition",  [](Bead * bed, pyarray pos) {bed->setPosition(to_vector(pos));})
        .def("radius", &Bead::radius)
        .def("radiusSqr", &Bead::radiusSqr)
        .def("resize", &Bead::resize)
        .def("volume", &Bead::volume)
        .def("dragCoefficient", &Bead::dragCoefficient)
        .def("next", &Bead::next, py::return_value_policy::reference)
        .def("prev", &Bead::prev, py::return_value_policy::reference)
        .def("toBead",  [](Object * obj) {return Bead::toBead(obj);},  py::return_value_policy::reference);
        
        
        
        
    
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
        .def("nbPoints", [](const Solid * sol) {return sol->nbPoints();})
        .def("toSolid",  [](Object * obj) {return Solid::toSolid(obj);},  py::return_value_policy::reference);
         
         /**
            @TODO : ADD SPECIALIZED SOLID CLASSES
         */
         
    py::class_<SolidProp,Property>(m, "SolidProp")
        .def_readwrite("drag", &SolidProp::drag)
        .def_readwrite("viscosity", &SolidProp::viscosity)
        .def_readwrite("steric", &SolidProp::steric)
        .def_readwrite("steric_range", &SolidProp::steric_range)
        .def_readwrite("confine", &SolidProp::confine)
        .def_readwrite("confine_stiffness", &SolidProp::confine_stiffness)
        .def_readwrite("confine_space", &SolidProp::confine_space)
        .def_readwrite("display", &SolidProp::display)
        .def_readwrite("display_fresh", &SolidProp::display_fresh);
        
}

