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
class simul;

/// A class containing a vect of fibers with the same properties
class FiberGroup : public std::vector<Fiber*>
{
    public:
        FiberProp * prop;
        FiberGroup() = default;
        FiberGroup(FiberProp * p) : FiberGroup() {prop = p ;};
        ~FiberGroup() = default;
};

/// A dictionary of FiberGroups
typedef std::map<std::string, FiberGroup> Fibers;

/// A time frame ; basicaly 
struct Frame 
{
    
    py::dict objects;
    Fibers fibers;
    int time;
};

//Frame & prepare_frame(int ) ;
Frame & prepared_frame( Simul * , int ) ;

#endif

/*
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
*/