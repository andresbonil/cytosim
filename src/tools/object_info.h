#ifndef PYTHON_OBJECT_H
#define PYTHON_OBJECT_H

#include <iostream>
#include <string>
#include <stack>
#include <map>
#include <unordered_map>
#include <vector>
#include "real.h"
#include "vector.h"





//typedef std::vector<ObjReport> obj_vec;
typedef std::unordered_map<std::string,real> prop_reals;

typedef std::unordered_map<std::string,std::string> prop_strings;
//typedef std::unordered_map<std::string,ObjReport> map_objs;

typedef std::unordered_map<std::string,std::string> string_dict;
typedef std::unordered_map<std::string,real> real_dict;
typedef std::unordered_map<std::string,int> int_dict;
typedef std::unordered_map<std::string,Vector> vector_dict;
typedef std::unordered_map<std::string,real_array> array_dict;

/// Basically a dictionnary
struct ObjectInfo {

    real_dict reals;
    int_dict ints;
    vector_dict vecs;
    string_dict strings;

};

#endif
