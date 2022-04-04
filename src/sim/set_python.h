#ifndef PYTHON_SET_H
#define PYTHON_SET_H

//#include "objecter.h"
#include "objecter_python.h"


typedef std::vector<ObjReport *> reportSet;

class Object;
class ObjectSet;

struct SetReport {
    
    real_dict reals;
    int_dict ints;
    vector_dict vecs;
    array_dict arrays;
    string_dict strings;
    reportSet reports;
    
    
    SetReport() ;
    ~SetReport()  = default;
};

/*
struct PPyObj : public PyObj {
    
    using PyObj::PyObj;
    
};
*/

//py::class_<PyObj>(m,"PyObj");

//py::class_<PyObj>(m, "PyObj").attr("id",py::cast(&PyObj::id));

#endif
