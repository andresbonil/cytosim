// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
/**
 This is a program to analyse simulation results:
 it reads a trajectory-file, and print some data from it.
*/


/**

  To use in python : move the cytosim...._.so file to a folder with *.cmo files
    We recommend using cym/aster.cym for a demo.
   
    Then run : 
        
    import cytosim
    sim = cytosim.open()
    sim.prop.timestep
    frame = cytosim.frame(0)
    fibers = frame["microtubule"]
    fibers.prop.segmentation = 1.337    # <- Yes, yes, yes.
    fibers[0].points()
    fibers[0].id()
    fibers[0].prop
    core = frame["core"][0]
    core.points()
     
    # etc...
*/

#include "report_python.h"

namespace py = pybind11;

//Simul simul;
FrameReader reader;
int status = -2;

//extern Simul simul;
extern FrameReader reader;
extern int status;


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
    //arg.set(verbose, "verbose");
    if ( arg.use_key("-") ) verbose = 0;

    try
    {
        RNG.seed();
        sim->loadProperties();
        reader.openFile(input);
        Cytosim::all_silent();
        status = -1;

    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        return nullptr;
    }

    return sim;
}


/**
 * @brief Prepares a frame from the simulation
 * @param sim
 * @param frame
 * @return 
 */
Frame * prepare_frame( Simul * sim, int frame) 
{
    reader.loadFrame(*sim, frame);
    Frame * current = new Frame;
    
    // First deal with the fibers
    PropertyList plist = sim->properties.find_all("fiber");
    for ( Property * i : plist )
    {
        FiberProp * fp = static_cast<FiberProp*>(i);
        current->fibers[fp->name()] = FiberGroup(fp);
    }
    for (auto fiber = sim->fibers.first(); fiber != sim->fibers.last() ; fiber = fiber->next() ) {
        current->fibers[fiber->property()->name()].push_back(fiber);
    }
    current->fibers[sim->fibers.last()->property()->name()].push_back(sim->fibers.last());
    
    // Then deals with spaces
    plist = sim->properties.find_all("space");
    for ( Property * i : plist )
    {
        SpaceProp * fp = static_cast<SpaceProp*>(i);
        current->spaces[fp->name()] = SpaceGroup(fp);
    }
    for (auto space = sim->spaces.first(); space != sim->spaces.last() ; space = space->next() ) {
        current->spaces[space->property()->name()].push_back(space);
    }
    current->spaces[sim->spaces.last()->property()->name()].push_back(sim->spaces.last());
    
    // Deals with solids
    plist = sim->properties.find_all("solid");
    for ( Property * i : plist )
    {
        SolidProp * fp = static_cast<SolidProp*>(i);
        current->solids[fp->name()] = SolidGroup(fp);
    }
    for (auto solid = sim->solids.first(); solid != sim->solids.last() ; solid = solid->next() ) {
        current->solids[solid->property()->name()].push_back(solid);
    }
    current->solids[sim->solids.last()->property()->name()].push_back(sim->solids.last());
    
    // Deal with whatever
    // ....
    
    // Grouping in dict
    for (const auto &[name, group] : current->fibers) {
        current->objects[py::cast(name)] = group;
    }
    for (const auto &[name, group] : current->spaces) {
        current->objects[py::cast(name)] = group;
    }
    for (const auto &[name, group] : current->solids) {
        current->objects[py::cast(name)] = group;
    }
    
    //return *current;
    return current;
}

int get_status() {
    return status;
}



/**
 * @brief  A module to get cytosim in python 
 * @return 
 * @TODO : a lot should be put in specific files
 */
PYBIND11_MODULE(cytosim, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    /// Loading properties into the module
    load_prop_classes(m);
        
    /// Python interface to Object
    py::class_<Object>(m, "Object")
        .def("id",  [](const Object * obj) {return obj->identity();})
        .def("info",  [](const Object * obj) {return to_dict(obj->info());});
    
    load_fiber_classes(m);
    load_solid_classes(m);    
    
    /// Python interface to FiberGroup : behaves roughly as a Python list of fibers
    py::class_<FiberGroup>(m, "FiberGroup")
        .def("__len__", [](const FiberGroup &v) { return v.size(); })
        .def_readwrite("prop",   &FiberGroup::prop , py::return_value_policy::reference)
        .def("__iter__", [](FiberGroup &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("__getitem__",[](const FiberGroup &v, size_t i) {
                 if (i >= v.size()) {
                     throw py::index_error();
                 }
                 return v[i];
             }, py::return_value_policy::reference);
             
    /// Python interface to SolidGroup : behaves roughly as a Python list of solids
    py::class_<SolidGroup>(m, "SolidGroup")
        .def("__len__", [](const SolidGroup &v) { return v.size(); })
        //.def("prop",  [](const SolidGroup & v) {return to_dict(v.prop->info());})
        .def("__iter__", [](SolidGroup &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("__getitem__",[](const SolidGroup &v, size_t i) {
                 if (i >= v.size()) {
                     throw py::index_error();
                 }
                 return v[i];
             }, py::return_value_policy::reference);
             
    /// Python interface to SpaceGroup, given as a single space
    py::class_<SpaceGroup>(m, "Space")
        //.def("__len__", [](const FiberGroup &v) { return v.size(); })
        //.def("prop",  [](const SpaceGroup & v) {return to_dict(v.prop->info());})
        .def("info",  [](const SpaceGroup & v) {return to_dict(v[0]->info());});     
        
    /// Python interface to timeframe : behaves roughly as a Python dict of ObjectGroup
    py::class_<Frame>(m, "Timeframe")
        .def_readwrite("fibers", &Frame::fibers, py::return_value_policy::reference)
        .def("__iter__", [](Frame &f) {
            return py::make_iterator(f.objects.begin(), f.objects.end());
        }, py::keep_alive<0, 1>())
        .def("__getitem__",[](const Frame &f, std::string s) {
                 return f.objects[py::cast(s)];
             }, py::return_value_policy::reference);
    
    /// Python interface to simul
    py::class_<Simul>(m, "Simul")
        .def_readwrite("prop",   &Simul::prop , py::return_value_policy::reference)
        .def("frame", [](Simul * sim, size_t i) 
            {return prepare_frame(sim, i);}, py::return_value_policy::reference);
        /*
       .def("__getitem__",[](Simul * sim, size_t i) {
                 //if (i >= v.size()) {
                 //    throw py::index_error();
                 //}
                 return prepared_frame(sim, i);
             });*/
    
    m.def("status", &get_status, "Status of the simul  : loaded or not");
    m.def("open", &open, "loads simulation from object files", py::return_value_policy::reference);
    //m.def("frame", &prepare_frame, "load frame");
    
}



//Frame & prepare_frame(int frame) {
//    return prepared_frame(&simul, frame);
//}


/*
 
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
 
  
  
 void select_frame(int frame) {
    reader.loadFrame(simul, frame);
    status = frame;
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
    //arg.set(verbose, "verbose");
    if ( arg.use_key("-") ) verbose = 0;

    try
    {
        RNG.seed();
        simul.loadProperties();
        reader.openFile(input);
        Cytosim::all_silent();
        status = -1;

    }
    catch( Exception & e )
    {
        std::clog << "Aborted: " << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return 0;
}


   
  */






    //m.def("get_first_fiber", &get_first_fiber, "load simulation");
    //m.def("select_frame", &select_frame, "load simulation");
    //m.def("report_frame_single", &report_fframe, "blaaa");
    //m.def("report_frame", &report_frame, "blaaa");
    //m.def("report_framer", &report_framer, "blaaa");
    ///m.def("report_fibers_frame", &report_framer2, "blaaa");
    //m.def("fibers", &get_fiber_set, "load simulation");
    //m.def("get_reals", &get_props, "A function that reports fiber frame f");

    /*
      auto a = py::class_<PyObj>(m, "object")
        .def_readwrite("id", &PyObj::id);
    a.def_readwrite("points", &PyObj::points);
    a.def_readwrite("props", &PyObj::props);
    
    auto b = py::class_<PySet>(m, "set");
    b.def_readwrite("objects", &PySet::objects);
    b.def_readwrite("props", &PySet::props);
    
    auto c = py::class_<PySetter>(m, "setter");
    //c.def_readwrite("props", &PySetter::props);
    
    auto d = py::class_<PyObjs>(m, "objects")
        .def_readwrite("props", &PyObjs::props)
        .def(py::init<>())
        .def("__len__", [](const PyObjs &v) { return v.size(); })
        .def("__iter__", [](PyObjs &v) {
            return py::make_iterator(v.begin(), v.end());
        }, py::keep_alive<0, 1>())
        .def("__getitem__",[](const PyObjs &v, size_t i) {
                 if (i >= v.size()) {
                     throw py::index_error();
                 }
                 return v[i];
             });
    py::class_<FiberSet>(m, "fiberSet")
        .def("__len__", [](const FiberSet * v) { return v->size(); })
        .def("__iter__", [](FiberSet *v) {
            return py::make_iterator(v->first(), v->last());
        }, py::keep_alive<0, 1>());
     */



/*

Fiber * get_first_fiber(int frame) {
    if (status >= -1 ) {
        reader.loadFrame(simul, frame);
        
        Fiber * fib = simul.fibers.first();
        return fib;
    }
}


PyObjs report_framer2(int frame) {
    if (status >= -1 ) {
        reader.loadFrame(simul, frame);
        SetReport * rep = simul.fibers.report();
        PyObjs fibers(rep);
        delete rep;
        return fibers;
    }
    else {
        return PyObjs();
    }
}


PyObj report_fframe(int fr) {
    if (status >= -1 ) {
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
    if (status >= -1 ) {
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
    if (status >= -1 ) {
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


PySet::PySet(SetReport* rep) {
        
        props = py::cast(rep->reals);
        props.attr("update")(py::cast(rep->strings));
        props.attr("update")(py::cast(rep->ints));
        props.attr("update")(py::cast(rep->vecs));
        
        for (auto obj: rep->objects) {
            objects.attr("append")(PyObj(obj));
        }
}

PySetter::PySetter(SetReport* rep) {
        
        props = py::cast(rep->reals);
        props.attr("update")(py::cast(rep->strings));
        props.attr("update")(py::cast(rep->ints));
        props.attr("update")(py::cast(rep->vecs));
        
        for (auto obj: rep->objects) {
            this->attr("append")(PyObj(obj));
        }
}

PyObjs::PyObjs(SetReport* rep) {
        
        props = py::cast(rep->reals);
        props.attr("update")(py::cast(rep->strings));
        props.attr("update")(py::cast(rep->ints));
        props.attr("update")(py::cast(rep->vecs));
        
        for (auto obj: rep->objects) {
            this->push_back(PyObj(obj));
        }
}



PyObj::PyObj(ObjReport* rep) {
                
        id = rep->id;
        
        points = pyarray(std::get<1>(rep->points),std::get<2>(rep->points), std::get<0>(rep->points));
        
        props = py::cast(rep->reals);
        props.attr("update")(py::cast(rep->strings));
        props.attr("update")(py::cast(rep->ints));
        props.attr("update")(py::cast(rep->vecs));
}
*/

/// returns an unordered map of PyObj
//map_objs get_objs(int i, int j) {
//    map_objs objs{{std::to_string(i),get_frame_obj(i)},{std::to_string(j),get_frame_obj(j)}};
 //   return objs;
//}