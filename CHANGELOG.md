# Pytosim 
## Compilation
First install pybind and then compile "report_python". Pytosim requires python >= 3.7 and a compiler supporting C++17.

```bash
$ python3 -m pip install -U --user pybind11
$ make -j4 report_python
```
This should yield a file cytosim.(...).so in your bin folder. E.g. : "cytosim.cpython-37m-x86_64-linux-gnu.so"

## Principle
Pytosim in an interface to native cytosim objects. However, there is currently still a sorting needed. Cytosim objects are better accessed through a Frame object. A frame is a dictionary of cytosim objects.

For example   

 ```python
    frame = sim.frame()
    print(frame.keys())
    fibers = frame["microtubule"]
```
Here fibers is a (python) list of (cytosim) Fiber objects. You can use native cytosim function for fibers, e.g.

```python
    fibers[0].nbPoints() 
```
Will yield the number of points.  
Additionally, a points() function has been defined :  
 
```python
     fibers[0].points()
```  
Yields a numpy array. 

## To load existinig sim:
Assuming that the cmo files and cytosim.-.so are in the current folder : 

```python
    import cytosim
    sim = cytosim.open()
    sim.prop.timestep 
    frame = cytosim.load(0)
    fibers = frame["microtubule"]
    fibers[0].points()
    fibers[0].id()
    fibers[0].join(fibers[1]) # <- yes, indeed
    core = frame["core"][0]
    core.points()
    while frame.loaded:
        print(frame.time)
        frame = frame.next()
```

## To run a simulation, from python
```python
    sim = cytosim.start('cym/aster.cym')
    frame = sim.frame() 
    fibers = frame['microtubule'] 
    fibers[0].join(fibers[1])    # <- Yes, yes, yes. 
    sim.step()
    sim.solve() 
```

# What changed 
Basically no code change was performed in cytosim except :   
- object.cc/h was changed to objecter.cc/h
    -> all files with "#include object.h" need to change to "#include objecter.h"  
- node.cc/h was changed to noder.cc/h  
     -> all files with "#include node.h" need to change to "#include noder.h"  
- In "sim_thread.cc", line 440 was commented : "//glApp::flashText0(str);"  
- makefile.inc and tools/makefile.inc were changed to allow compilation.  
- Then a lot of files were added to /tools  


