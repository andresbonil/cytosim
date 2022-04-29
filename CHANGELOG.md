# Cythosim
Cythosim is a C-Python build of cytosim.
## Compilation
First install pybind and then compile "report_python". Cythosim requires python >= 3.7 and a compiler supporting C++17.

```bash
$ python3 -m pip install -U --user pybind11
$ make -j4 report_python
```
This should yield a file cytosim.(...).so in your bin folder. E.g. : "cytosim.cpython-37m-x86_64-linux-gnu.so"

## Principle
Cythosim is an interface to native cytosim objects.

 ```python
    import cytosim
    sim = cytosim.start("cym.aster.cym")
    for fiber in sim.fibers:
        print(fiber.length())
```
Here sim.fibers is a (cytosim) FiberSet, i.e. the set of all fibers.

For complex simulations, with different types of fibers, couples, etc.,  we also offer a Frame object, in which objects are sorted by name : a Frame is a (python) dictionary of lists of cytosim objects, that all have the same property.

For example   

 ```python
    frame = sim.frame()
    mts = frame["microtubule"]
    for fiber in mts:
        print(fiber.length())
```
Here mts is a (python) list of (cytosim) Fiber objects. You can use native cytosim function for fibers, e.g.

```python
    mt = mts[0]
    mt.nbPoints()
```
Will yield the number of points.  
Additionally, a points() function has been defined, yielding a numpy array :  

```python
     mt.points()
```  

To know the methods available from an object, type dir():

```python
    print(dir(sim))
    print(dir(frame))
    print(dir(frame["core"][0]))
```

We can easily change property :
```python
    sim.prop.time_step = 0.1
    sim.prop.complete(sim)
    print(sim.prop.time_step)
    mts.prop.change_str("rigidity = 0.1", sim)
    print(mts.prop.rigidity)
```

## To load existing sim:
Assuming that the cmo files and cytosim.-.so are in the current folder :

```python
    import cytosim
    sim = cytosim.open()
    sim.prop.timestep
    frame = cytosim.load(0)
    fibers = frame["microtubule"]
    fibers[0].points()
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
