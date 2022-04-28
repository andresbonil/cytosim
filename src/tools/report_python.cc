// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and provides a python interface to it.
*/


/**

  To use in python : move the cytosim...._.so file to a folder with *.cmo files
    We recommend using cym/aster.cym for a demo.
   
    Then run : 
        
    import cytosim
    sim = cytosim.open()
    sim.prop.timestep 
    frame = cytosim.load(0)
    fibers = frame["microtubule"]
    fibers.prop.segmentation = 1.337    # <- Yes, yes, yes.
    fibers[0].points()
    fibers[0].id()
    fibers[0].join(fibers[1]) # <- yes, indeed
    core = frame["core"][0]
    core.points()
    while frame.loaded:
        print(frame.time)
        frame = frame.next()
     
    # etc...
*/

/*
@TODO : an interface for FiberSet (problem : cannot iterate because of FiberSet interface)
@TODO : support input arguments
 */
#include "report_python.h"

namespace py = pybind11;

/// Using global vars, sorry not sorry.
FrameReader reader;
bool __is_loaded__ = 0;
//SimThread * thread;
extern FrameReader reader;
extern bool __is_loaded__;
//extern thread;

void bar(void) {};

/// Open the simulation from the .cmo files
Simul * open()
{   
    
    int verbose = 1;
    int prefix = 0;
    
    Glossary arg;

    std::string input = TRAJECTORY;
    std::string str;

    Simul * sim = new Simul;
    
    unsigned period = 1;

    arg.set(input, ".cmo") || arg.set(input, "input");    
    if ( arg.use_key("-") ) verbose = 0;

    try
    {
        RNG.seed();
        sim->loadProperties();
        reader.openFile(input);
        Cytosim::all_silent();
        __is_loaded__ = 1;
    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        return nullptr;
    }

    return sim;
}

Simul * start(std::string fname) {
    std::cout << fname << std::endl;
    int n = fname.length();
    char inp[n] ;
    std::strcpy(inp, fname.c_str());
    Glossary arg;
    arg.read_string(inp,2);
    
    if ( ! arg.use_key("+") )
    {
        Cytosim::out.open("messages.cmo");
        Cytosim::log.redirect(Cytosim::out);
        Cytosim::warn.redirect(Cytosim::out);
    }
    
    Simul * simul = new Simul;
    std::cout << "new simul" << std::endl;
    try {
        simul->initialize(arg);
        std::cout << "initialized" << std::endl;
    }
    catch( Exception & e ) {
        print_magenta(std::cerr, e.brief());
        std::cerr << '\n' << e.info() << '\n';
    }
    catch(...) {
        print_red(std::cerr, "Error: an unknown exception occurred during initialization\n");
    }
    
    std::cout << "before parse" << std::endl;
    //arg.print_warning(std::cerr, 1, " on command line\n");
    time_t sec = TicToc::seconds_since_1970();
    //Parser(*simul, 1, 1, 1, 1, 1).readConfig();
    std::string file = simul->prop->config_file;
    std::string setup = file;
    std::cout << "before parse" << std::endl;
    Parser(*simul, 0, 1, 0, 0, 0).readConfig();
    //void foo(void) {};
    void (*foofoo)(void) = &bar;
    std::cout << "before thread" << std::endl;
    SimThread * thread = new SimThread(*simul, foofoo);
    std::cout << "made thread" << std::endl;
    //thread->period(1);
    thread->start();
    __is_loaded__ = 2;
    return simul;
}

/**
 * @brief Prepares a frame from the simulation
 * @param sim
 * @param frame
 * @return 
 */


/**
 * @brief  A module to get cytosim in python 
 * @return 
 * @TODO : a lot should be put in specific files
 */
PYBIND11_MODULE(cytosim, m) {
    m.doc() = "sim = cytosim.open() \n"
                "sim.prop.timestep \n"
                "frame = cytosim.load(0) \n"
                "fibers = frame['microtubule'] \n"
                "fibers.prop.segmentation = 1.337    # <- Yes, yes, yes. \n"
                "fibers[0].points() \n"
                "fibers[0].id() \n"
                "core = frame['core'][0] \n"
                "core.points() \n"
                "while frame.loaded: \n"
                "    print(frame.time) \n"
                "    frame = frame.next()"; // optional module docstring
        
    /// Loading properties into the module
    load_object_classes(m);
    auto pysim = load_simul_classes(m);
    load_fiber_classes(m);
    load_hand_classes(m);
    load_solid_classes(m);
    load_space_classes(m);
    load_single_classes(m);
    load_couple_classes(m);
    load_organizer_classes(m);
    
    /// We declare object groups
    // We can later add additional def to any of these groups
    auto fibs = declare_group(m, ObjGroup<Fiber,FiberProp>(), "FiberGroup");
    auto sols = declare_group(m, ObjGroup<Solid,SolidProp>(), "SolidGroup");
    auto spas = declare_group(m, ObjGroup<Space,SpaceProp>(), "SpaceGroup");
    auto beds = declare_group(m, ObjGroup<Bead,BeadProp>(), "BeadGroup");
    auto sfrs = declare_group(m, ObjGroup<Sphere,SphereProp>(), "SphereGroup");
    auto orgs = declare_group(m, ObjGroup<Organizer,Property>(), "OrganizerGroup");
    auto sins = declare_group(m, ObjGroup<Single,SingleProp>(), "SingleGroup");
    auto cous = declare_group(m, ObjGroup<Couple,CoupleProp>(), "CoupleGroup");
    
    /// Python interface to timeframe : behaves roughly as a Python dict of ObjectGroup
    py::class_<Frame>(m, "Timeframe")
        .def_readwrite("fibers", &Frame::fibers, py::return_value_policy::reference)
        .def_readwrite("time", &Frame::time)
        .def_readwrite("index", &Frame::index)
        .def_readwrite("loaded", &Frame::loaded)
        .def("next", [](Frame &f) {return prepare_frame(f.simul, &reader, f.index+1);}) //py::return_value_policy::reference
        .def("__iter__", [](Frame &f) {
            return py::make_iterator(f.objects.begin(), f.objects.end());
        }, py::keep_alive<0, 1>())
        .def("keys", [](Frame &f) {  return f.objects.attr("keys")() ; })
        .def("items", [](Frame &f) { return f.objects.attr("items")() ; })
        .def("__getitem__",[](const Frame &f, std::string s) {
                 return f.objects[py::cast(s)];
             }, py::return_value_policy::reference);

            
    /// Python interface to simul
    pysim.def("load", [](Simul * sim, size_t i) 
            {if (__is_loaded__) {return prepare_frame(sim, &reader, i);} ; return new Frame; }); //py::return_value_policy::reference
    pysim.def("frame", [](Simul * sim) 
            {if (__is_loaded__) {return make_frame(sim);} ; return new Frame; }); //py::return_value_policy::reference
            
    /// Opens the simulation from *.cmo files
    m.def("open", &open, "loads simulation from object files", py::return_value_policy::reference);
    m.def("start", &start, "loads simulation from config files", py::return_value_policy::reference);
    
}

