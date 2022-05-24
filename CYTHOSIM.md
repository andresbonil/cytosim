# Cythosim
Cythosim is an *experimental* C-Python interface for cytosim.

It allows to run, save, load simulations, as well as interact with native cytosim objects.   
Please refer to iPython notebooks for examples. The functions available for cytosim objects mirror the C++ functions and thus are detailed in cytosim documentation. Currently all functions are not yet implemented thou.
## Compilation
First install pybind and then compile "cythosim". Cythosim requires python >= 3.7 and a compiler supporting C++17. Code below should yield a file  
 cytosim.(...).so in your bin folder. E.g. : "cytosim.cpython-37m-x86_64-linux-gnu.so"

```bash
$ python3 -m pip install -U --user pybind11
$ make -j4 cythosim
```

## Principle
Cythosim is an interface to native cytosim objects.
### Running simulations

 ```python
    import cytosim
    sim = cytosim.start("cym.aster.cym")
    for fiber in sim.fibers:
        print(fiber.length())
    sim.run(10) # runs 10 timesteps
    sim.save()
```

### Loading simulations
Assuming that the cmo files and cytosim.-.so are in the current folder :

```python
    import cytosim
    sim = cytosim.open()
    while sim.next():
	time = sim.time()
	n_fibs = len(sim.fibers)
        print(f"At time {time}s, simul has {n_fibs} fibers")
```

### Simulation frame
For complex simulations, with different types of fibers, couples, etc.,  we also offer a Frame object, in which objects are sorted by name : a Frame is a (python) dictionary of lists of cytosim objects, that all have the same property. For example  :

 ```python
    frame = sim.frame()
    mts = frame["microtubule"]
    print(len(mts))
    frame = frame.next()
```
While it takes extra time to create a frame, it may be worth it.


## What changed
Very little code change was performed in cytosim except :   
- node.cc/h was changed to noder.cc/h  
     -> all files with "#include node.h" need to change to "#include noder.h"  
- In "sim_thread.cc", line 440 was commented : "//glApp::flashText0(str);"   
- makefile.inc and tools/makefile.inc were changed to allow compilation.   
- Then a lot of files were added to /tools  and cpython/  
- Glossary can now export mTerms through the public function Glossary::terms()  
- In simul.h : added Simul::prepare_meca : a wrapper for sMeca.prepare  
- In simul.h and simul_solve.cc : addition of Simul::prepared_solve() : basically Simul::solve without sMeca.prepare.  
- In mecable.h : addition of mecable::nonConstData() : like data() but not const.  

## How to use Cythosim on other branches of cytosim ?
If you want to use cythosim on other branches of cytosim, you can "easily" do so with a few operations :  
- Download cythosim in an other folder.  
- Copy paste files src/tools/cythosim.* and src/tools/makefile.inc into your own src/tools.  
- Copy folder src/cpython into your own src/  
- (possibly) Copy makefile.inc into your own folder.  
- Rename "src/base/node.cc/h" to "src/base/noder.cc/h" and change all occurences of "#include 'node.h'" to "#include 'noder.h'" in your code.  
- Copy "src/base/glossary.h" to your own "src/base". What is important is line 226 of cythosim's glossary.h.  
- Copy "src/sim/simul.h" and "src/sim/simul_solve.cc" to your own src/sim. What is important is the declaration of Simul::prepare_meca and Simul::prepared_solve  
- Copy "src/sim/mecable.h" to your own "src/sim/" folder. What is important is the declaration of nonConstData().  
