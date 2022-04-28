# Pytosim !
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

Or, in live mode ! 

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
- node.* was changed to noder.*  
     -> all files with "#include node.h" need to change to "#include noder.h"  
- In "sim_thread.cc", line 440 was commented : "//glApp::flashText0(str);"  
- makefile.inc and tools/makefile.inc were changed to allow compilation.  
- Then a lot of files were added to /tools  


