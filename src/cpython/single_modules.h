#include "single.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Single;
class Object;
/// Converts an object to a single if possible
static Single* toSingle(Object * obj)
{
    if ( obj  &&  obj->tag() == 's' )
        return static_cast<Single*>(obj);
    return nullptr;
}

/// a utility to enrich the cytosim python module
void load_single_classes(py::module_ &m) {
     /// Python interface to single
    py::class_<Single,Object>(m, "Single")
        .def("toSingle",  [](Object * s) {return toSingle(s);}, py::return_value_policy::reference)
        .def("state", &Single::state)
        .def("__getitem__",[](Single *s, int i) { // We can call Single[0]  to get the first hand ! thus couple[0].attachEnd(...) is available
            if (i==0) {return s->hand();} else {throw py::index_error(); }
            return (Hand*)nullptr;} 
            , py::return_value_policy::reference)
        .def("__len__",  [](Single * s) {return (int)1;})
        .def("hand", &Single::hand, py::return_value_policy::reference)
        .def("attached", &Single::attached)
        .def("fiber", &Single::fiber, py::return_value_policy::reference)
        .def("abscissa", &Single::abscissa)
        .def("posHand",  [](Single * s) {return to_numpy(s->posHand());})
        .def("dirFiber",  [](Single * s) {return to_numpy(s->dirFiber());})
        .def("attach", &Single::attach)
        .def("attachEnd",  [](Single * s, Fiber *f, int end) {return s->attachEnd(f,static_cast<FiberEnd>(end));})
        .def("moveToEnd",  [](Single * s, int end) {return s->moveToEnd(static_cast<FiberEnd>(end));})
        .def("detach", &Single::detach)
        .def("position",  [](Single * s) {return to_numpy(s->position());})
        .def("mobile", &Single::mobile)
        .def("translate",  [](Single * s, pyarray vec) 
            {   Vector p = to_vector(vec);
                return s->translate(p);})
        .def("setPosition",  [](Single * s, pyarray vec) 
            {   Vector p = to_vector(vec);
                return s->setPosition(p);})
        .def("randomizePosition", &Single::randomizePosition)
        .def("posFoot",  [](Single * s) {return to_numpy(s->posFoot());})
        .def("sidePos",  [](Single * s) {return to_numpy(s->sidePos());})
        .def("base", &Single::base)
        .def("mobile", &Single::mobile)
        .def("force",  [](Single * s) {return to_numpy(s->force());})
        .def("next", &Single::next)
        .def("prev", &Single::prev)
        .def("confineSpace", &Single::confineSpace);
        
        

    py::class_<SingleProp,Property>(m, "SingleProp")
        .def_readwrite("hand", &SingleProp::hand)
        .def_readwrite("stiffness", &SingleProp::stiffness)
        .def_readwrite("length", &SingleProp::length)
        .def_readwrite("diffusion", &SingleProp::diffusion)
        .def_readwrite("fast_diffusion", &SingleProp::fast_diffusion)
        .def_readwrite("confine", &SingleProp::confine)
        .def_readwrite("confine_space", &SingleProp::confine_space)
        .def_readwrite("activity", &SingleProp::activity)
        .def_readwrite("hand_prop", &SingleProp::hand_prop);
    
     py::class_<SingleSet,ObjectSet>(m, "SingleSet")
		.def("__getitem__",[](SingleSet * set, int i) {
				int s = set->size();
                if (i<0) {i+=s;} // Python time negative indexing
				if (i >= s or i<0) {
					 throw py::index_error();
				}
				Single * obj = set->firstID();
				while (i) {
					--i; // I know this is slow, but ...
					obj = set->nextID(obj); 
				}
				return obj;
             }, py::return_value_policy::reference);
}

