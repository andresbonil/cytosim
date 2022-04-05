// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/
#include "report_python.h"

namespace py = pybind11;

Simul simul;
FrameReader reader;
int status = 0;

extern Simul simul;
extern FrameReader reader;
extern int status;

PyObj::PyObj(ObjReport* rep) {
                
        id = rep->id;
        
        points = py::array_t<double>(std::get<1>(*rep->points),std::get<2>(*rep->points), std::get<0>(*rep->points));
        
        props = py::cast(*rep->reals);
        props.attr("update")(py::cast(*rep->strings));
        props.attr("update")(py::cast(*rep->ints));
        props.attr("update")(py::cast(*rep->vecs));        
}

PySet::PySet(SetReport* rep) {
        
        props = py::cast(*rep->reals);
        props.attr("update")(py::cast(*rep->strings));
        props.attr("update")(py::cast(*rep->ints));
        props.attr("update")(py::cast(*rep->vecs));
        
        for (auto obj: *rep->objects) {
            objects.attr("append")(PyObj(obj));
        }
}

PySetter::PySetter(SetReport* rep) {
        
        props = py::cast(*rep->reals);
        props.attr("update")(py::cast(*rep->strings));
        props.attr("update")(py::cast(*rep->ints));
        props.attr("update")(py::cast(*rep->vecs));
        
        for (auto obj: *rep->objects) {
            this->attr("append")(PyObj(obj));
        }
}



int load_simul()
{   
    int verbose = 1;
    int prefix = 0;
    
    Glossary arg;

    std::string input = TRAJECTORY;
    std::string str;

    
    unsigned period = 1;

    arg.set(input, ".cmo") || arg.set(input, "input");
    arg.set(verbose, "verbose");
    if ( arg.use_key("-") ) verbose = 0;

    try
    {
        RNG.seed();
        simul.loadProperties();
        reader.openFile(input);
        Cytosim::all_silent();
        status = 1;

    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return 0;
}

PyObj report_fframe(int fr) {
    if (status == 1 ) {
        reader.loadFrame(simul, fr);
        ObjReport * rep = simul.fibers.firstID()->report();
        PyObj obj(rep);
        delete rep;
        return obj;
    }
    else {
        return PyObj();
        }
}

PySet report_frame(int frame) {
    if (status == 1 ) {
        reader.loadFrame(simul, frame);
        SetReport * rep = simul.fibers.report();
        PySet fibers(rep);
        delete rep;
        return fibers;
    }
    else {
        return PySet();
    }
}

PySetter report_framer(int frame) {
    if (status == 1 ) {
        reader.loadFrame(simul, frame);
        SetReport * rep = simul.fibers.report();
        PySetter fibers(rep);
        delete rep;
        return fibers;
    }
    else {
        return PySetter();
    }
}



/// A function that reads a frame and returns a numpy array 
// Returns the first position of the first fiber
pyarray report_loaded_frame(int fr)
{
    

    unsigned frame = fr;
    //std::array<real,DIM> pts;

    if (status == 1 ) {
        reader.loadFrame(simul, frame);
        
        //const real * data = simul.fibers.firstID()->data();
        //int size = simul.fibers.firstID()->nbPoints();
        
        
       
        //Vector pos = simul.fibers.firstID()->posP(0);

        //for (unsigned p=0;p<DIM;++p) {
        //      pts[p] = pos[p];
        //    }
        //std::vector<int> sizes = {size, (int)DIM};
        //std::vector<int> strides = {DIM*sizeof(real), sizeof(real)};
        //real_array array_pts{data, sizes, strides};
        //pyarray test = pyarray ({std::get<1>(array_pts),DIM},{DIM*8,8},std::get<0>(array_pts),py::return_value_policy::copy);
        //py::array_t<real> test = py::array_t<double>(std::get<1>(array_pts),std::get<2>(array_pts), std::get<0>(array_pts));
    
        ObjReport * rep = simul.fibers.firstID()->report();
        
        int id = rep->id;
        
        
        py::array_t<real> test = py::array_t<double>(std::get<1>(*rep->points),std::get<2>(*rep->points), std::get<0>(*rep->points));
        return test;
    }
    else {
        std::array<real,1> pts{0};
        pyarray pypts;
        return pypts;
    }


    // Reading reporting the first position of the first fiber
    //pyarray pypts = py::cast(pts);
    
}


int get_status() {
    return status;
}



/// Showcasing making dictionaries
py::dict get_props() {
    //prop_reals reals{std::pair<std::string,real>{"TEST",1.0}};
    prop_reals reals{{"TEST",1.0},{"TEST2",2.0}};
    prop_strings strings{{"type","essai"}};
    py::dict dict;
    dict = py::cast(reals);
    dict.attr("update")(py::cast(strings));
    return dict;
}

/// returns an unordered map of PyObj
//map_objs get_objs(int i, int j) {
//    map_objs objs{{std::to_string(i),get_frame_obj(i)},{std::to_string(j),get_frame_obj(j)}};
 //   return objs;
//}

PYBIND11_MODULE(cytosim, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    //auto a = py::class_<ObjReport>(m, "ObjReport")
    //    .def_readwrite("id", &ObjReport::id);
        //.def_readwrite("pos0", &PyObj::position);
    //a.def_readwrite("pos0", &ObjReport::pos0);
    auto a = py::class_<PyObj>(m, "object")
        .def_readwrite("id", &PyObj::id);
    a.def_readwrite("points", &PyObj::points);
    a.def_readwrite("props", &PyObj::props);
    
    auto b = py::class_<PySet>(m, "set");
    b.def_readwrite("objects", &PySet::objects);
    b.def_readwrite("props", &PySet::props);
    
    auto c = py::class_<PySetter>(m, "setter");
    
    //c.def_readwrite("props", &PySetter::props);
    
    m.def("get_reals", &get_props, "A function that reports fiber frame f");
    m.def("status", &get_status, "blaaa");
    m.def("report_loaded", &report_loaded_frame, "blaaa");
    m.def("report_frame_single", &report_fframe, "blaaa");
    m.def("report_frame", &report_frame, "blaaa");
    m.def("report_framer", &report_framer, "blaaa");
    m.def("load", &load_simul, "load simulation");
}


/*

  To use in python : move the example...._.so file to a folder with *.cmo 
   
    Then run : 
    import example
    example.add(1,2)
    example.getobj(1)
    etc...
*/
