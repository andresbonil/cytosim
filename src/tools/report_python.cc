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
int verbose = 1;
int prefix = 0;
size_t cnt = 0;

typedef py::array_t<double> pyarray;
typedef std::vector<PyObj> obj_vec;
typedef std::unordered_map<std::string,real> prop_reals;
typedef std::unordered_map<std::string,std::string> prop_strings;
typedef std::unordered_map<std::string,PyObj> map_objs;

/// A function that reads a frame and returns a numpy array 
// Returns the first position of the first fiber
pyarray report_fiber_frame(int fr)
{
    
    Glossary arg;

    std::string input = TRAJECTORY;
    std::string str;

    unsigned frame = fr;
    unsigned period = 1;

    arg.set(input, ".cmo") || arg.set(input, "input");
    arg.set(verbose, "verbose");
    if ( arg.use_key("-") ) verbose = 0;

    Simul simul;
    FrameReader reader;

    try
    {
        RNG.seed();
        simul.loadProperties();
        reader.openFile(input);
        Cytosim::all_silent();

        reader.loadFrame(simul, frame);
    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        //return EXIT_FAILURE;
    }


    // Reading reporting the first position of the first fiber
    std::array<real,DIM> pts;
    Vector pos = simul.fibers.firstID()->posP(0);
    for (unsigned p=0;p<DIM;++p) {
        pts[p] = pos[p];
    }
    pyarray pypts = py::cast(pts);


    /// check if all specified parameters were used:
    arg.print_warning(std::cerr, cnt, "\n");
    
    return pypts;
}

// Just returns a PyObj, i.e. a struct with a py::array and a real
PyObj report_fib_frame(int fr)
{
    
    Glossary arg;

    std::string input = TRAJECTORY;
    std::string str;

    unsigned frame = fr;
    unsigned period = 1;

    arg.set(input, ".cmo") || arg.set(input, "input");
    arg.set(verbose, "verbose");
    if ( arg.use_key("-") ) verbose = 0;

    Simul simul;
    FrameReader reader;

    try
    {
        RNG.seed();
        simul.loadProperties();
        reader.openFile(input);
        Cytosim::all_silent();

        reader.loadFrame(simul, frame);
    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        //return EXIT_FAILURE;
    }

    PyObj pobj(simul.fibers.firstID());

    return pobj;
}

// Simplest example
int add(int i, int j) {
    return i + j;
}

// reports a py::array directy
pyarray ff_report(int i) {
    // Returns a numpy array for the position of the first point of the fiber at frame i
    pyarray arr = report_fiber_frame(i);
    return arr;
}

// returns a py_obj
PyObj get_frame_obj(int i) {
    PyObj fib = report_fib_frame(i);
    return fib;
}

// returns a vector of PyObj
obj_vec get_frames(int i, int j) {
    obj_vec fibs{get_frame_obj(i), get_frame_obj(j) };
    return fibs;
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
map_objs get_objs(int i, int j) {
    map_objs objs{{std::to_string(i),get_frame_obj(i)},{std::to_string(j),get_frame_obj(j)}};
    return objs;
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    auto a = py::class_<PyObj>(m, "PyObj")
        .def_readwrite("id", &PyObj::id);
        //.def_readwrite("pos0", &PyObj::position);
    a.def_readwrite("pos0", &PyObj::position);

    m.def("add", &add, "A function that adds two numbers");
    m.def("report", &ff_report, "A function that reports fiber frame f");
    m.def("getobj", &get_frame_obj, "A function that reports fiber frame f");
    m.def("get_two", &get_frames, "A function that reports fiber frame f");
    m.def("get_reals", &get_props, "A function that reports fiber frame f");
    m.def("get_objs", &get_objs, "A function that reports fiber frame f");
}
