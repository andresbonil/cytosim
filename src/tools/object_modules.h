#include "objecter.h"
#include <pybind11/pybind11.h>
#include "python_utilities.h"
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;

class Object;

/// a utility to enrich the cytosim python module
void load_object_classes(py::module_ &m) {
	py::class_<ObjectSet>(m, "ObjectSet")
		.def("add",  [](ObjectSet * set, Object * obj) {return set->add(obj) ;})
		.def("remove",  [](ObjectSet * set, Object * obj) {return set->remove(obj) ;})
		.def("link",  [](ObjectSet * set, Object * obj) {return set->link(obj) ;})
		.def("unlink",  [](ObjectSet * set, Object * obj) {return set->unlink(obj) ;})
		.def("erase",  [](ObjectSet * set, Object * obj) {return set->erase(obj) ;})
		.def("erase_all",  [](ObjectSet * set) {return set->erase() ;})
		.def("size",  [](ObjectSet * set) {return set->size() ;})
		.def("__len__",  [](ObjectSet * set) {return set->size() ;})
		.def("first",  [](ObjectSet * set) {return set->first() ;}, py::return_value_policy::reference)
		.def("last",  [](ObjectSet * set) {return set->last() ;}, py::return_value_policy::reference)
		.def("findID",  [](ObjectSet * set, int n) 
			{return set->findID((ObjectID) n) ;}, py::return_value_policy::reference)
		.def("findObject",  [](ObjectSet * set, Property * p)
			{return set->findObject(p) ;}, py::return_value_policy::reference)
		.def("__getitem__",[](ObjectSet * set, int i) {
				int s = set->size();
				if (i >= s) {
					 throw py::index_error();
				}
				Object * obj = set->first();
				
				while (i) {
					--i;
					obj = obj->next();
				}
				return obj;
             }, py::return_value_policy::reference);	
	
     /// Python interface to Organizer
    py::class_<Object>(m, "Object")
        .def("reference",  [](Object * obj) {return obj->reference() ;})
        .def("property",  [](Object * obj) {return obj->property() ;})
        .def("position", [](const Object * obj) {return to_numpy(obj->position());})
        .def("next",  [](Object * obj) {return obj->next() ;}, py::return_value_policy::reference)
		.def("__next__",  [](Object * obj) {return obj->next() ;}, py::return_value_policy::reference)
        .def("prev",  [](Object * obj) {return obj->prev() ;}, py::return_value_policy::reference)
        .def("id",  [](const Object * obj) {return obj->identity();})
        .def("points", [](const Object * obj) {return pyarray();});
        
}
