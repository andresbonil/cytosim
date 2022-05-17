#include "fiber.h"
#include "classic_fiber.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Fiber;
class Object;

/// a utility to enrich the cytosim python module
void load_fiber_classes(py::module_ &m) {
    /// Python interface to Mecable
    py::class_<Chain,Mecable>(m, "Chain")
        .def("nbSegments",  [](Chain * chn) {return chn->nbSegments();})
        .def("lastSegment",  [](Chain * chn) {return chn->lastSegment();})
        .def("setStraight",  [](Chain * chn, pyarray pos, pyarray dir) 
            {   Vector p = to_vector(pos);
                Vector d = to_vector(dir);
                return chn->setStraight(p,d);})
        .def("setStraightLength",  [](Chain * chn, pyarray pos, pyarray dir, real len) 
            {   Vector p = to_vector(pos);
                Vector d = to_vector(dir);
                return chn->setStraight(p,d,len);})
        .def("placeEnd",  [](Chain * chn, int end) {return chn->placeEnd((FiberEnd)end);})
        .def("setEquilibrated",  [](Chain * chn, real len, real persil) {return chn->setEquilibrated(len, persil);})
        .def("birthTime",  [](Chain * chn) {return chn->birthTime();})
        .def("age",  [](Chain * chn) {return chn->age();})
        .def("exactEnd",  [](Chain * chn, int a) {
            return new Mecapoint(chn->exactEnd((FiberEnd)a)); })
        .def("interpolateEndM", [](Chain * chn) {
            return new Interpolation(chn->interpolateEndM());})
        .def("interpolateEndP",  [](Chain * chn) {
            return new Interpolation(chn->interpolateEndP()); })
        .def("interpolateCenter", [](Chain * chn) {
            return new Interpolation(chn->interpolateCenter()); })
        .def("interpolateEnd",  [](Chain * chn, int a)   {
            return new Interpolation(chn->interpolateEnd((FiberEnd)a)); })
        .def("interpolateFromEnd",  [](Chain * chn, real ab, int a) {
            return new Interpolation(chn->interpolate(ab,(FiberEnd)a)); })
        .def("interpolate",  [](Chain * chn, real ab) {
            return new Interpolation(chn->interpolate(ab)); })
        .def("length",  [](Chain * chn) {return chn->length();})
        .def("length1",  [](Chain * chn) {return chn->length1();})
        .def("trueLength", &Chain::trueLength)
        .def("betweenMP",  [](Chain * chn, real a) {return chn->betweenMP(a);})
        .def("outsideMP",  [](Chain * chn, real a) {return chn->outsideMP(a);})
        .def("belowP",  [](Chain * chn, real a) {return chn->belowP(a);})
        .def("aboveM",  [](Chain * chn, real a) {return chn->aboveM(a);})
        .def("whichEndDomain",  [](Chain * chn, real a, real b) {return chn->whichEndDomain(a,b);})
        .def("clampedIndexM",  [](Chain * chn, real a) {return chn->clampedIndexM(a);})
        .def("setOrigin",  [](Chain * chn, real a) {return chn->setOrigin(a);})
        .def("abscissaM",  [](Chain * chn) {return chn->abscissaM();})
        .def("abscissaC",  [](Chain * chn) {return chn->abscissaC();})
        .def("abscissaP",  [](Chain * chn) {return chn->abscissaP();})
        .def("abscissaPoint",  [](Chain * chn, real a) {return chn->abscissaPoint(a);})
        .def("abscissaEnd",  [](Chain * chn, int a) {return chn->abscissaPoint((FiberEnd)a);})
        .def("abscissaFrom",  [](Chain * chn, real dis, int a) {return chn->abscissaFrom(dis,(FiberEnd)a);})
        .def("someAbscissa",  [](Chain * chn, real dis, int a, int mod, real alf) {return chn->someAbscissa(dis,(FiberEnd)a, mod, alf);})
        .def("posM",  [](Chain * chn, real a) {return to_numpy(chn->posM(a));})
        .def("pos",  [](Chain * chn, real a) {return to_numpy(chn->pos(a));})
        .def("posFrom",  [](Chain * chn, real a, int ref) {return chn->posFrom(a,(FiberEnd)ref);})
        .def("posMiddle",  [](Chain * chn) {return to_numpy(chn->posMiddle());})
        .def("posEnd",  [](Chain * chn, int end) {return to_numpy(chn->posEnd((FiberEnd)end));})
        .def("posEndP",  [](Chain * chn) {return to_numpy(chn->posEndP());})
        .def("posEndM",  [](Chain * chn) {return to_numpy(chn->posEndM());})
        .def("netForceEndM",  [](Chain * chn) {return to_numpy(chn->netForceEndM());})
        .def("netForceEndP",  [](Chain * chn) {return to_numpy(chn->netForceEndP());})
        .def("projectedForceEndM",  [](Chain * chn) {return chn->projectedForceEndM();})
        .def("projectedForceEndP",  [](Chain * chn) {return chn->projectedForceEndP();})
        .def("projectedForceEndM",  [](Chain * chn, int end) {return chn->projectedForceEnd((FiberEnd) end);})
        .def("avgDirection",  [](Chain * chn) {return to_numpy(chn->avgDirection());})
        .def("segmentation",  [](Chain * chn) {return chn->segmentation();})
        .def("flipChainPolarity",  [](Chain * chn) {return chn->flipChainPolarity();})
        .def("curvature",  [](Chain * chn, unsigned p) {return chn->curvature(p);})
        .def("bendingEnergy0",  [](Chain * chn) {return chn->bendingEnergy0();})
        .def("planarIntersect",  [](Chain * chn, unsigned s, pyarray vec, real a) {return chn->planarIntersect(s, to_vector(vec), a);})
        .def("growM",  [](Chain * chn, real a) {return chn->growM(a);})
        .def("addSegmentM",  [](Chain * chn) {return chn->addSegmentM();})
        .def("cutM",  [](Chain * chn, real a) {return chn->cutM(a);})
        .def("growP",  [](Chain * chn, real a) {return chn->growP(a);})
        .def("addSegmentP",  [](Chain * chn) {return chn->addSegmentP();})
        .def("cutP",  [](Chain * chn, real a) {return chn->cutP(a);})
        .def("grow",  [](Chain * chn, int ref, real a) {return chn->grow((FiberEnd)ref,a);})
        .def("adjustLength",  [](Chain * chn, real a, int ref) {return chn->adjustLength(a,(FiberEnd)ref);})
        .def("truncateM",  [](Chain * chn, unsigned a) {return chn->truncateM(a);})
        .def("truncateP",  [](Chain * chn, unsigned a) {return chn->truncateP(a);});        
    
    py::class_<Mecafil,Chain>(m, "Mecafil")
        .def("tension",  [](Mecafil * mec, unsigned p) {return mec->tension(p);})
        .def("dragCoefficient",  [](Mecafil * mec) {return mec->dragCoefficient();})
        .def("leftoverMobility",  [](Mecafil * mec) {return mec->leftoverMobility();});
    
    py::class_<Fiber,Mecafil>(m, "Fiber")
        .def("points",  [](Fiber * fib) {return get_obj_points(fib);})
        .def("toFiber",  [](Object * obj) {return Fiber::toFiber(obj);},  py::return_value_policy::reference)
        .def("nbPoints",  [](Fiber * fib) {return fib->nbPoints();})
        .def("cutM",  [](Fiber * fib, real len) {return fib->cutM(len);})
        .def("cutP",  [](Fiber * fib, real len) {return fib->cutP(len);})
        .def("sever",  [](Fiber * fib, real a, int p, int m) {return fib->sever(a, p, m);})
        .def("severNow",  [](Fiber * fib, real a) {return fib->severNow(a);})
        .def("severP",  [](Fiber * fib, real a) {return fib->severP(a);})
        .def("nowSever",  [](Fiber * fib) {return fib->severNow();})
        .def("join",  [](Fiber * fib, Fiber * fob) {return fib->join(fob);})
        .def("stateM",  [](const Fiber * fib) {return py::cast(fib->dynamicStateM());})
        .def("stateP",  [](const Fiber * fib) {return py::cast(fib->dynamicStateP());})
        .def("setStateM",  [](Fiber * fib, int stat) {return fib->setDynamicStateM(stat);})
        .def("setStateP",  [](Fiber * fib, int stat) {return fib->setDynamicStateP(stat);})
        .def("setDynamicState",  [](Fiber * fib, int end, int stat) {return fib->setDynamicState(static_cast<FiberEnd>(end),stat);})
        .def("freshAssembly",  [](Fiber * fib, int end) {return fib->freshAssembly(static_cast<FiberEnd>(end));})
        .def("nbHands",  [](Fiber * fib) {return fib->nbHands();})
        .def("nbHandsInRange",  [](Fiber * fib, real amin, real amax, int end) {return fib->nbHandsInRange(amin, amax, static_cast<FiberEnd>(end));})
        .def("__next__", [](const Fiber * fib) {return fib->next();});
        
        
    py::class_<ClassicFiber,Fiber>(m, "ClassicFiber")
        .def("freshAssemblyM",  [](ClassicFiber * fib) {return fib->freshAssemblyM();})
        .def("freshAssemblyP",  [](ClassicFiber * fib) {return fib->freshAssemblyP();});
        /**
         * 
            @TODO : complete with fiber base functions
         * 
        */
         
         /**
            @TODO : ADD SPECIALIZED FIBER CLASSES
         */
         
    /// Python interface to FiberProp
    py::class_<FiberProp,Property>(m, "FiberProp")
        .def_readwrite("segmentation", &FiberProp::segmentation)
        .def_readwrite("rigidity", &FiberProp::rigidity)
        .def_readwrite("min_length", &FiberProp::min_length)
        .def_readwrite("max_length", &FiberProp::max_length)
        .def_readwrite("total_polymer", &FiberProp::total_polymer)
        .def_readwrite("persistent", &FiberProp::persistent)
        .def_readwrite("viscosity", &FiberProp::viscosity)
        .def_readwrite("drag_radius", &FiberProp::drag_radius)
        .def_readwrite("drag_length", &FiberProp::drag_length)
        .def_readwrite("drag_model", &FiberProp::drag_model)
        .def_readwrite("drag_gap", &FiberProp::drag_gap)
        .def_readwrite("binding_key", &FiberProp::binding_key)
        .def_readwrite("lattice", &FiberProp::lattice)
        .def_readwrite("lattice_unit", &FiberProp::lattice_unit)
        .def_readwrite("confine", &FiberProp::confine)
        .def_readwrite("confine_space", &FiberProp::confine_space)
        .def_readwrite("confine_stiffness", &FiberProp::confine_stiffness)
        .def_readwrite("steric", &FiberProp::steric)
        .def_readwrite("steric_radius", &FiberProp::steric_radius)
        .def_readwrite("steric_range", &FiberProp::steric_range)
        .def_readwrite("glue", &FiberProp::glue)
        .def_readwrite("glue_single", &FiberProp::glue_single)
        .def_readwrite("activity", &FiberProp::activity)
        .def_readwrite("display_fresh", &FiberProp::display_fresh)
        .def_readwrite("display", &FiberProp::display);
        
        py::enum_<FiberEnd>(m,"FiberEnd")
            .value("NO_END", NO_END)
            .value("PLUS_END", PLUS_END)
            .value("MINUS_END", MINUS_END)
            .value("BOTH_ENDS", BOTH_ENDS)
            .value("ORIGIN", ORIGIN)
            .value("CENTER", CENTER)
            .export_values();

        py::enum_<AssemblyState>(m,"AssemblyState")
            .value("STATE_WHITE", STATE_WHITE)
            .value("STATE_GREEN", STATE_GREEN)
            .value("STATE_YELLOW", STATE_YELLOW)
            .value("STATE_ORANGE", STATE_ORANGE)
            .value("STATE_RED", STATE_RED)
            .value("STATE_BLACK", STATE_BLACK)
            .export_values();
        
}

