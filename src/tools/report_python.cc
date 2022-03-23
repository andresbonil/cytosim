// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/

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

size_t cnt = 0;
Simul simul;
FrameReader reader;

extern Simul simul;
extern FrameReader reader;

extern int status;
int status = 0;
typedef py::array_t<double> pyarray;
typedef std::vector<PyObj> obj_vec;
typedef std::unordered_map<std::string,real> prop_reals;
typedef std::unordered_map<std::string,std::string> prop_strings;
typedef std::unordered_map<std::string,PyObj> map_objs;

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


/// A function that reads a frame and returns a numpy array 
// Returns the first position of the first fiber
pyarray report_loaded_frame(int fr)
{
    

    unsigned frame = fr;
    std::array<real,DIM> pts;

    if (status == 1 ) {
        reader.loadFrame(simul, frame);
        Vector pos = simul.fibers.firstID()->posP(0);

        for (unsigned p=0;p<DIM;++p) {
            pts[p] = pos[p];
        }
    }
    else {
        std::array<real,1> pts{0};
    }


    // Reading reporting the first position of the first fiber
    pyarray pypts = py::cast(pts);


    return pypts;
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

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    auto a = py::class_<PyObj>(m, "PyObj")
        .def_readwrite("id", &PyObj::id);
        //.def_readwrite("pos0", &PyObj::position);
    a.def_readwrite("pos0", &PyObj::position);

    
    m.def("get_reals", &get_props, "A function that reports fiber frame f");
    m.def("status", &get_status, "blaaa");
    m.def("report_loaded", &report_loaded_frame, "blaaa");
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