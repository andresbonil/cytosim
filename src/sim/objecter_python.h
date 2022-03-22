#ifndef PYTHON_OBJECT_H
#define PYTHON_OBJECT_H

#include "objecter.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;
class Object;

struct PyObj {
  Vector pos0;  
  py::array_t<double> position;
  int id;
  PyObj(const Object* obj  ) ;
  int get_id() {return id;} ;
 ~PyObj()  = default;
};

/*
struct PPyObj : public PyObj {
    
    using PyObj::PyObj;
    
};
*/

//py::class_<PyObj>(m,"PyObj");

//py::class_<PyObj>(m, "PyObj").attr("id",py::cast(&PyObj::id));

#endif