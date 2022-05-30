// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and provides a python interface to it.
   
   * @TODO :    - manage to return ObjectSet from simul, in order to not necessitate frame()
                - bead and sphere
                - live player + python ? o_O
                - specialized classes, including dynamic spaces
                - and so on and so forth
*/


/**

  To use in python : move the cytosim...._.so file to a folder with *.cmo files
    We recommend using cym/aster.cym for a demo.
   
 
    import cytoplay
    sim = cytoplay.start('cym/aster.cym')
    def runtimeCheck(simul):
        return simul.time()
    cytoplay.setRuntimeCheck(runtimeCheck)
    cytoplay.play()

*/

/*
@TODO : an interface for FiberSet (problem : cannot iterate because of FiberSet interface)
@TODO : support input arguments
 */
#include "cythosim.h"

#include "opengl.h"
#include "player.h"
#include "view.h"
#include "gle.h"
#include <pybind11/functional.h>
#include <thread>
namespace py = pybind11;

Player player;




/// Using global vars, sorry not sorry.
FrameReader reader ;
Simul&      simul = player.simul;
int __is_loaded__ = 0;
SimThread & thread = player.thread;
bool __saved__ = 0;
Parser * parser;
PlayerProp&  prop = player.prop;
DisplayProp& disp = player.disp;


#  include "glut.h"
#  include "glapp.h"
#  include "fiber_prop.h"
#  include "fiber_disp.h"
#  include "point_disp.h"
using glApp::flashText;
#  include "play_keys.cc"
#  include "play_menus.cc"
#  include "play_mouse.cc"


extern FrameReader reader;
extern int __is_loaded__;
extern SimThread & thread;
extern bool __saved__;
extern Parser * parser;
extern Player player;
extern PlayerProp& prop;
extern DisplayProp& disp;

/// A holder for normalKey callback
inline std::function<unsigned char(unsigned char, int, int)>& normalKey()
{
    // returns a different object for each thread that calls it
    static thread_local std::function<unsigned char(unsigned char, int, int)> fn;
    return fn;
}
/// A proxy for the normalKeyy callback
inline void proxyNormalKey(unsigned char c, int i, int j){ c = normalKey()(c, i ,j ); processNormalKey(c,i,j); };

inline std::function<int(int, int, const Vector3&, int)>& mouseClick()
{
    // returns a different object for each thread that calls it
    static thread_local std::function<int(int, int, const Vector3&, int)> mc;
    return mc;
}
/// A proxy for the normalKeyy callback
inline void proxyMouseClick(int i, int j, const Vector3& v, int k){int c = mouseClick()(i ,j, v, k );
    processMouseClick(i,j,v,c); };

/// A holder for runtime callback
inline std::function<void(Simul&)>& runtimeCheck()
{
    // returns a different object for each thread that calls it
    static thread_local std::function<void(Simul&)> rt;
    return rt;
}

/// Displays the simulation live
void displayLive(View& view)
{
    // Also adds a callback to an external function through caller->runtimeCheck
    if ( 0 == thread.trylock() )
    {
        // read and execute commands from incoming pipe:
        thread.readInput(32);
        //thread.debug("display locked");
        if ( simul.prop->display_fresh )
        {
            player.readDisplayString(view, simul.prop->display);
            simul.prop->display_fresh = false;
        }
        
        player.prepareDisplay(view, 1);
        player.displayCytosim();
        
        // external callback
        runtimeCheck()(simul);
        thread.unlock();
        
    }
    else
    {
        thread.debug("display: trylock failed");
        glutPostRedisplay();
    }
}

/// Initiates a simulation
Simul * start(std::string fname ) {
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
    
    try {
        simul.initialize(arg);
    }
    catch( Exception & e ) {
        print_magenta(std::cerr, e.brief());
        std::cerr << '\n' << e.info() << '\n';
    }
    catch(...) {
        print_red(std::cerr, "Error: an unknown exception occurred during initialization\n");
    }
    
    //arg.print_warning(std::cerr, 1, " on command line\n");
    time_t sec = TicToc::seconds_since_1970();
    
    std::string file = simul.prop->config_file;
    std::string setup = file;
    
    parser = new Parser(simul, 0, 1, 0, 0, 0);
    parser->readConfig();
    
    thread.period(prop.period);
    thread.start();
    __is_loaded__ = 2;

    // Default null callbacks
    normalKey() = [](unsigned char c, int i, int j) {return c;} ;
    mouseClick() = [](int i, int j, const Vector3 v, int k) {return k;} ;
    runtimeCheck() = [](Simul& sim) {};
    
    return &simul;
}

void play_default(std::string opt){
//#ifdef __APPLE__
#if (1)
    int argc = 1;
    std::string st = "";
    char * inp[1];
    std::strcpy(*inp, st.c_str());
    glutInit(&argc, inp);
#endif
    Glossary arg = Glossary(opt);
    
    glApp::setDimensionality(DIM);
    if ( arg.use_key("fullscreen") )
        glApp::setFullScreen(1);
    View& view = glApp::views[0];
    view.read(arg);
    disp.read(arg);
    simul.prop->read(arg);
    view.setDisplayFunc(displayLive);
    
    // Definining the callbacks
    glApp::actionFunc(proxyMouseClick);
    glApp::actionFunc(processMouseDrag);
    glApp::normalKeyFunc(proxyNormalKey);
    glApp::createWindow(displayLive);
    
    try
    {
        gle::initialize();
        player.setStyle(disp.style);
        rebuildMenus();
        glutAttachMenu(GLUT_RIGHT_BUTTON);
        glutMenuStatusFunc(menuCallback);
        if ( glApp::isFullScreen() )
            glutFullScreen();
        glutTimerFunc(200, timerCallback, 0);
    }
    catch ( Exception & e )
    {
        print_magenta(std::cerr, e.brief());
        std::cerr << '\n' << e.info() << '\n';
    }
    
    try
    {
        glutMainLoop();
    }
    catch ( Exception & e )
    {
        print_magenta(std::cerr, e.brief());
        std::cerr << '\n' << e.info() << '\n';
    }
}


/// A python module to run or play cytosim
PYBIND11_MODULE(cytoplay, m) {
    m.doc() = "# live mode only \n"
                "sim = cytoplay.start('cym/aster.cym') \n"
                "def runtimeCheck(simul): \n"
                "   print(simul.time()) \n"
                "cytoplay.setRuntimeCheck(runtimeCheck) \n"
                "cytoplay.play() \n";
                 // optional module docstring
    
    /// Loading properties into the module
    load_object_classes(m);
    load_meca_classes(m);
    load_point_classes(m);
    auto pysim = load_simul_classes(m);
    load_glossary_classes(m);
    load_solid_classes(m);
    load_fiber_classes(m);
    load_hand_classes(m);
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
        .def("update", [](Frame &f) {  return make_frame(f.simul) ; })
        .def("keys", [](Frame &f) {  return f.objects.attr("keys")() ; })
        .def("items", [](Frame &f) { return f.objects.attr("items")() ; })
        .def("__getitem__",[](const Frame &f, std::string s) {
                 return f.objects[py::cast(s)];
             }, py::return_value_policy::reference);

    /// Python interface to play/start a simulation
    m.def("simul",[](){return &simul ;}, py::return_value_policy::reference);
    m.def("start",[](std::string fname ){return start(fname) ;}, py::return_value_policy::reference);
    m.def("play", [](py::args args) {
        int nargs = args.size();
        if (nargs == 0) { play_default("")  ; }
        else {
            std::string opt;
            for (auto arg : args) {
                opt += py::cast<std::string>(arg);
                }
            std::cout << opt << std::endl;
            play_default(opt);
            }
        }, py::call_guard<py::gil_scoped_release>());
    m.def("setNormalKey",[](py::function f) {
        normalKey() = py::cast<std::function<unsigned char(unsigned char, int, int)>>(f);
        });
    m.def("setRuntimeCheck",[](py::function f) {
        runtimeCheck() = py::cast<std::function<void(Simul&)>>(f);
    });
    m.def("setMouseClick",[](py::function f) {
        mouseClick() = py::cast<std::function<int(int, int, Vector3, int)>>(f);
    });
 
    m.def("str_to_glos", &str_to_glos, "converts string to Glossary");
    

    /// Expading Python interface to simul
    pysim.def("frame", [](Simul * sim) 
            {return make_frame(sim);} ); //py::return_value_policy::reference
    pysim.def("writeObjects",  [](Simul * sim) {
        sim->writeObjects(sim->prop->trajectory_file,__saved__,1);
        if (!__saved__) {__saved__ = 1;}  ;} );
    pysim.def("writeObjects",  [](Simul * sim, std::string & str) {
        sim->writeObjects(str,__saved__,1);
        if (!__saved__) {__saved__ = 1;}  ;} );
    pysim.def("save", [](Simul * sim) {
            sim->writeObjects(sim->prop->trajectory_file,__saved__,1);
            if (!__saved__) {__saved__ = 1;};
            sim->writeProperties(&sim->prop->property_file[0],1);
        });
    pysim.def("add",  [](Simul * sim, py::args args) {
        std::string name = "";
        std::string how = "";
        int many = 1;
        int nargs = args.size();
        if (nargs>0) {
            name = py::cast<std::string>(args[0]);
        }
        if (nargs>1) {
            how = py::cast<std::string>(args[1]);
        }
        if (nargs>2) {
            many = py::cast<int>(args[2]);
        }
        Glossary glos = Glossary(how);
        ObjectList objects;
        for (int i=0;i<many;++i) {
            auto objs = parser->execute_new(name, glos);
            objects.push_back(objs[0]);
        }
        return objects;
        }, py::return_value_policy::reference);
    pysim.def("cut",  [](Simul * sim, std::string & name, std::string & where) {
            Glossary glos = Glossary(where);
            parser->execute_cut(name, glos);
            });
    pysim.def("delete",  [](Simul * sim, std::string & name, std::string & how, int number) {
            Glossary glos = Glossary(how);
            parser->execute_delete(name, glos, number);
            });
    pysim.def("import",  [](Simul * sim, std::string & file, std::string & what, std::string & how) {
            Glossary glos = Glossary(how);
            parser->execute_import(file, what, glos);
            });
    pysim.def("export",  [](Simul * sim, std::string & file, std::string & what, std::string & how) {
            Glossary glos = Glossary(how);
            parser->execute_export(file, what, glos);
            });
    pysim.def("run",  [](Simul * sim, py::args args) {
            std::string how = "";
            int many = 1;
            int nargs = args.size();
            if (nargs>0) {
                many = py::cast<int>(args[0]);
            }
            if (nargs>1) {
                how = py::cast<std::string>(args[1]);
            }
            Glossary glos = Glossary(how);
            parser->execute_run(many, glos,0);
            });
    pysim.def("set",  [](Simul * sim, std::string & cat, std::string & name, std::string & how) {
            Glossary glos = Glossary(how);
            parser->execute_set(cat, name, glos);
            return glos;
            });
}

