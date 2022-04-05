#include "objecter_python.h"

ObjReport::ObjReport() {
    // Object identity
    /* const Object* obj 
    id = obj->identity();
    ints.insert({"id",id});
    
    const real * data = nullptr;
    int size = 0;
    
    std::vector<int> sizes = {size, (int)DIM};
    std::vector<int> strides = {DIM*sizeof(real), sizeof(real)};
    points = real_array{data, sizes, strides};
     */
    id = -1;
    points = new real_array;
    
    reals  = new real_dict ;
    ints = new int_dict;
    vecs = new vector_dict;
    strings = new string_dict;

} ;

