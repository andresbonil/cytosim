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

// contains an adress and 
typedef std::tuple<real[], int> real_array;
typedef std::vector<PyObj> obj_vec;
typedef std::unordered_map<std::string,real> prop_reals;
typedef std::unordered_map<std::string,std::string> prop_strings;
typedef std::unordered_map<std::string,PyObj> map_objs;

typedef std::unordered_map<std::string,std::string> string_dict;
typedef std::unordered_map<std::string,real> real_dict;
typedef std::unordered_map<std::string,Vector> vector_dict;
typedef std::unordered_map<std::string,pyarray> pyarray_dict;
typedef std::unordered_map<std::string,real_array> array_dict;

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



int load_simul();
pyarray report_loaded_frame(int);
py::dict get_props();

#endif