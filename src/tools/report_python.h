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

typedef py::array_t<double> pyarray;


/// Reporter is a construction to report several kind of things to a python dictionarry
struct Reporter {
    vector_dict* vectors;
    array_dict* arrays;
    real_dict* reals;
    string_dict* strings;
    Reporter() {
        vectors = nullptr;
        arrays = nullptr;
        reals = nullptr;
        strings = nullptr;
    }
    ~Reporter() = default;
};

struct PyObj {
    pyarray points;
    int id;
    py::dict props;
};



int load_simul();
pyarray report_loaded_frame(int);
py::dict get_props();

#endif
