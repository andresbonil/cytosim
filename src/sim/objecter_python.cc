#include "objecter_python.h"

PyObj::PyObj(const Object* obj ) {
    pos0 = Vector(1,0);
    std::array<real,DIM> pts;
    for (int i=0;i<DIM;++i) {
        pts[i] = pos0[i];
    }
    //std::cout << pts << std::endl;
    position = py::cast(pts);
    id = obj->identity();
} ;

/*  
PyObj::~PyObj() {
    id = 0;
    pos0 = 0;
}
 */
