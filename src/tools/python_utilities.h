#ifndef UTILITIES_H
#define UTILITIES_H
#include "fiber.h"
#include <pybind11/pybind11.h>
//#include <pybind11/numpy.h>
//#include <pybind11/stl.h>
namespace py = pybind11;


typedef py::array_t<real> pyarray;

/// Converts a real array * to numpy array
pyarray & to_numpy(real_array * rar) {
    if (rar) {
        pyarray * arr =new pyarray(std::get<1>(*rar),std::get<2>(*rar), std::get<0>(*rar));
        return *arr;
    } else {
        pyarray * arr = new pyarray();
        return *arr;
    }
}

pyarray & to_numpy(Vector vec) {    
    pyarray * par = new pyarray;
#if DIM==1
    *par = py::cast(std::vector<real>{vec[0]});
#elif DIM==2
    *par = py::cast(std::vector<real>{vec[0],vec[1]});
#else
    *par = py::cast(std::vector<real>{vec[0],vec[1],vec[2]});
#endif
    return *par;
}

/// Converts an ObjectInfo * to a python dict
py::dict & to_dict(ObjectInfo * info) {
    py::dict * dico = new py::dict;
    dico->attr("update")(py::cast(info->reals));
    dico->attr("update")(py::cast(info->strings));
    dico->attr("update")(py::cast(info->ints));
    dico->attr("update")(py::cast(info->vecs));
    
    return *dico;
}


#endif