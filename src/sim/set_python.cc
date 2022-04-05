#include "set_python.h"

SetReport::SetReport() {
    //for ( Object * obj=Set->first(); obj; obj=obj->next() )
    //{
    //    reports.push_back(obj->report());
    //}
     
    
    reals  = new real_dict ;
    ints = new int_dict;
    vecs = new vector_dict;
    strings = new string_dict;
    objects = new reportSet;
}

/*  
PyObj::~PyObj() {
    id = 0;
    pos0 = 0;
}
 */
