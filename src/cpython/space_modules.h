#include "space.h"
#include "common.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Space;
class Object;
/// Converts an object to a space if possible;
static Space* toSpace(Object * obj)
{
    if ( obj  &&  obj->tag() == 'e' )
        return static_cast<Space*>(obj);
    return nullptr;
}

/// a utility to enrich the cytosim python module
void load_space_classes(py::module_ &m) {
     /// Python interface to space
    py::class_<Space,Object>(m, "Space")
        .def("thickness", &Space::thickness)
        .def("resize",  [](Space * sp, std::string sizes) {
            Glossary glos = Glossary(sizes); sp->resize(glos) ;})
        .def("volume", &Space::volume)
        .def("boundaries",  [](const Space * sp)
            {Vector V,W; 
            sp->boundaries(V,W);
                return std::vector<pyarray>{to_numpy(V),to_numpy(W)};
                })
        .def("inside", [](const Space * sp, pyarray pos) {return sp->inside(to_vector(pos));})
        .def("project", [](const Space * sp, pyarray pos) {return to_numpy(sp->project(to_vector(pos)));})
        .def("allInside", [](const Space * sp, pyarray pos, real rad) {return sp->allInside(to_vector(pos),rad);})
        .def("allOutside", [](const Space * sp, pyarray pos, real rad) {return sp->allOutside(to_vector(pos),rad);})
        .def("toSpace",  [](Object * s) {return toSpace(s);}, py::return_value_policy::reference)
        .def("max_extension", &Space::max_extension)
        .def("projectDeflated", [](const Space * sp, pyarray pos, real rad) 
            {return to_numpy(sp->projectDeflated(to_vector(pos),rad));})
        .def("distanceToEdgeSqr", [](const Space * sp, pyarray pos) {return sp->distanceToEdgeSqr(to_vector(pos));})
        .def("distanceToEdge", [](const Space * sp, pyarray pos) {return sp->distanceToEdge(to_vector(pos));})
        .def("signedDistanceToEdge", [](const Space * sp, pyarray pos) {return sp->signedDistanceToEdge(to_vector(pos));})
        .def("bounce", [](const Space * sp, pyarray pos) {return to_numpy(sp->bounce(to_vector(pos)));})
        .def("normalToEdge", [](const Space * sp, pyarray pos) {return to_numpy(sp->normalToEdge(to_vector(pos)));})
        .def("randomPlace", [](const Space * sp) {return to_numpy(sp->randomPlace());})
        .def("randomPlaceNearEdge", [](const Space * sp, real rad, int tries) {return to_numpy(sp->randomPlaceNearEdge(rad,tries));})
        .def("randomPlaceOnEdge", [](const Space * sp, real rad, int tries) {return to_numpy(sp->randomPlaceNearEdge(rad,tries));})
        .def("estimateVolume", &Space::estimateVolume);

    py::class_<SpaceProp,Property>(m, "SpaceProp");
            
    py::enum_<Confinement>(m,"Confinement")
        .value("CONFINE_OFF", CONFINE_OFF)
        .value("CONFINE_INSIDE", CONFINE_INSIDE)
        .value("CONFINE_OUTSIDE", CONFINE_OUTSIDE)
        .value("CONFINE_ON", CONFINE_ON)
        .value("CONFINE_ALL_INSIDE", CONFINE_ALL_INSIDE)
        .value("CONFINE_ALL_OUTSIDE", CONFINE_ALL_OUTSIDE)
        .value("CONFINE_PLUS_END", CONFINE_PLUS_END)
        .value("CONFINE_MINUS_END", CONFINE_MINUS_END)
        .value("CONFINE_BOTH_ENDS", CONFINE_BOTH_ENDS)
        .value("CONFINE_PLUS_OUT", CONFINE_PLUS_OUT)
        .value("CONFINE_POINT", CONFINE_POINT)
        .value("CONFINE_RANGE", CONFINE_RANGE)
        .export_values();

}

