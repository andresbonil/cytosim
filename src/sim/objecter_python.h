#ifndef PYTHON_OBJECT_H
#define PYTHON_OBJECT_H

#include <iostream>
#include <string>
#include <stack>
#include <map>
#include <unordered_map>
#include <vector>
#include "objecter.h"



typedef std::vector<int> int_vect;
// contains adress, sizes, and strides
typedef std::tuple<const real*, int_vect, int_vect> real_array;
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

//class Object;
//class ObjectSet;
/*
struct ObjReport {
    int id;
    
    real_array points;
    
    real_dict reals;
    int_dict ints;
    vector_dict vecs;
    string_dict strings;
    
  
  
  //ObjReport(const Object* obj  ) ;
  ObjReport() = default;
  ~ObjReport()  = default;
};
*/
#endif
