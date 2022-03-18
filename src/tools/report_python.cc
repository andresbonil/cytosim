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

namespace py = pybind11;
int verbose = 1;
int prefix = 0;
size_t cnt = 0;

typedef py::array_t<double> pyarray;

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



int add(int i, int j) {
    return i + j;
}

pyarray ff_report(int i) {
    // Returns a numpy array for the position of the first point of the fiber at frame i
    pyarray arr = report_fiber_frame(i);
    return arr;
}

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function that adds two numbers");
    m.def("report", &ff_report, "A function that reports fiber frame f");
}
