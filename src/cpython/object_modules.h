#include "object.h"
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
                if (i<0) {i+=s;} // Python time negative indexing
				if (i >= s or i<0) {
					 throw py::index_error();
				}
				Object * obj = set->first();
				while (i) {
					--i; // I know this is slow, but ...
					obj = obj->next(); // Maybe objectSet should derive from std::vect ?
				}
				return obj;
             }, py::return_value_policy::reference)
        .def("addObject",  [](ObjectSet * set, py::args args) {
            std::string name = "";
            std::string how = "";
            int many = 1;
            //int nargs = py::cast<int>args.attr("__len__")();
            int nargs = args.size();
            if (nargs>0) {
                name = py::cast<std::string>(args[0]);
            }
            if (nargs>1) {
                how = py::cast<std::string>(args[1]);
            }
            if (nargs>2) {
                many = py::cast<int>(args[2]);
            }
            Glossary glos = Glossary(how);
            ObjectList objects;
            for (int i=0;i<many;++i) {
                auto objs = set->newObjects(name, glos);
                set->simul.add(objs);
                objects.push_back(objs[0]);
            }
            return objects;
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
    
    
    py::class_<ObjectList>(m, "ObjectList")
        .def("__getitem__",[](ObjectList & lis, int i) {
				int s = lis.size();
                if (i<0) {i+=s;} // Python time negative indexing
				if (i >= s or i<0) {
					 throw py::index_error();
				}
				Object * obj = lis[i];
				return obj;
             }, py::return_value_policy::reference)
        .def("__len__", [](const ObjectList &v) { return v.size(); })
        /*.def("__getitem__",[](ObjectList & lis, int i) {
				int s = lis.size();
				if (i >= s) {
					 throw py::index_error();
				}
				Object * obj = lis[i];
				return obj;
             }, py::return_value_policy::reference)
         .def("__iter__", [](ObjectList & lis) {
            return py::make_iterator(lis.begin(), lis.end());
        }, py::keep_alive<0, 1>())
        .def("size", &ObjectList::size)
        .def("capacity", &ObjectList::capacity)
        .def("empty", &ObjectList::empty);*/
        
}
