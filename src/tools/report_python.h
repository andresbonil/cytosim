#ifndef REPORT_PYTHON_H
#define REPORT_PYTHON_H

#include <fstream>
#include <sstream>

#include "stream_func.h"
#include "frame_reader.h"
#include "iowrapper.h"
#include "glossary.h"
#include "messages.h"
#include "splash.h"
#include "parser.h"
#include "simul.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "objecter_python.h"
#include "objecter_python.h"
namespace py = pybind11;

typedef py::array_t<real> pyarray;

/// Reporter is a construction to report several kind of things to a python dictionarry


/// PyObj is the python representation of an object report
struct PyObj {
    pyarray points;
    int id;
    py::dict props;
    PyObj(ObjReport* );
    PyObj() = default;
    ~PyObj() = default;
};

/// PySet is the python representation of an set report
struct PySet {
    
    py::dict props;
    py::list objects;
    PySet(SetReport* );
    PySet() = default;
    ~PySet() = default;
};

/// PySetter is the python representation of a set report
class PySetter : public py::list
{
    public:
        py::dict props;
        PySetter() = default;
        PySetter(SetReport*);
        ~PySetter() = default;
};

/// PyObjs is the python representation of a set report 
class PyObjs : public std::vector<PyObj>
{
    public:
        py::dict props;
        PyObjs() = default;
        PyObjs(SetReport*);
        ~PyObjs() = default;
};



#endif
