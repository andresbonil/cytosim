// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef SIM_THREAD_H
#define SIM_THREAD_H

#include <pthread.h>
#include "simul.h"
#include "parser.h"
#include "frame_reader.h"


/// SimThread is used to run a simulation in a dedicated thread
class SimThread : private Parser
{
    /// disabled default constructor
    SimThread();
    
private:
    
    /// Simulation object
    Simul           simul;

    /// reader used to access frames in a trajectory file
    FrameReader     reader;

    
    /// callback invoked when the thread is halted
    void           (*hold_callback)(void);
    
    /// slave thread
    pthread_t       mChild;
    
    /// a flag reflecting if child thread is running or not
    bool            mState;
    
    /// a flag to indicate that child thread should terminate or restart
    int             mFlag;
    
    /// mutex protecting write access to simulation state
    pthread_mutex_t mMutex;

    /// condition variable used to control the thread execution
    pthread_cond_t  mCondition;
    
    /// counter for hold()
    unsigned int    mHold;
    
    /// period for hold()
    unsigned int    mPeriod;

    
    /// the current Single being controlled with the mouse
    mutable Single * mHandle;
    
    /// return the SingleProp used for the handles
    SingleProp *  getHandleProperty() const;

    /// make a new SingleProp for the handles with given attachment range
    SingleProp *  makeHandleProperty(real range);
    
    /// return list of Handles
    ObjectList    allHandles(SingleProp const*) const;

    /// True if thread is 'mChild'
    bool          isChild() const { return pthread_equal(pthread_self(), mChild); }
    
public:
    
    /// run the simulation live
    void          run();
    
    /// continue to run a simulation beyond its normal termination
    void          extend_run();

    /// redefines Interface::hold(), will be called repeatedly during parsing
    void          hold();
    
    /// print message to identify slave
    void          debug(const char *) const;
    
    /// create a SimThread with given holding function callback
    SimThread(void (*callback)(void));
    
    /// destructor
    ~SimThread();

#if ( 1 )

    /// lock access to the Simulation data
    void       lock()    { pthread_mutex_lock(&mMutex); }
    
    /// unlock access to the Simulation data
    void       unlock()  { pthread_mutex_unlock(&mMutex);}
    
    /// try to lock access to the Simulation data
    int        trylock() { return pthread_mutex_trylock(&mMutex); }

    /// unlock access to data and wait for the condition
    int        wait()    { return pthread_cond_wait(&mCondition, &mMutex); }
    
    /// send signal to other threads
    void       signal()  { if ( mState ) pthread_cond_signal(&mCondition); }

#else
    
    /// lock access to the Simulation data
    void       lock()    { pthread_mutex_lock(&mMutex); debug("lock"); }
    
    /// unlock access to the Simulation data
    void       unlock()  { pthread_mutex_unlock(&mMutex); debug("unlock"); }
    
    /// try to lock access to the Simulation data
    int        trylock() { int R=pthread_mutex_trylock(&mMutex); debug(R?"trylock failed":"trylock"); return R; }
    
    /// wait for the condition to be
    int        wait()    { debug(" wait"); return pthread_cond_wait(&mCondition, &mMutex); }
    
    /// signal other thread to continue
    void       signal()  { debug(" signal"); pthread_cond_signal(&mCondition); }
    
#endif
    
    /// Simul reference
    Simul&     sim() { return simul; }
    
    /// set how many 'hold()' are necessary to halt the thread
    void       period(unsigned int c) { mPeriod = c; }
    
    /// true if child thread is running
    bool       alive() const { return mState; }
    
    /// terminate (slave) thread
    void       terminate();
    
    /// start the thread that will run a simulation
    void       start();
    
    /// continue to run the simulation after its normal termination
    int        extend();
    
    /// perform one simulation step
    void       step();
    
    /// gently stop the simulation
    void       stop();

    /// stop the simulation
    void       cancel();
    
    /// restart engine
    void       restart();

    /// clear the simulation world
    void       clear();
    
    /// execute commands from standard input, return number of lines processed
    size_t     readInput(size_t max_nb_lines);
    
    /// halt the live simulation, read the config file and change the object parameters
    void       reloadParameters(std::string const& file);
    
    /// execute given code
    int        execute(std::string const&, std::string const&);
    
    /// export simulation Propertes and Objects to file
    void       exportObjects(bool binary);
    
    /// export properties to file
    void       writeProperties(std::ostream&, bool prune);

    
    /// open trajectory file for input
    void       openFile(std::string const& name) { reader.openFile(name); }
    
    /// true if ready to read from file
    bool       goodFile()  const { return reader.good(); }
    
    /// status of file
    int        eof()       const { return reader.eof(); }
    
    /// rewind file
    void       rewind()          { lock(); reader.rewind(); unlock(); }
 
    /// index of current frame
    int        currFrame() const { return reader.currFrame(); }
    
    /// attempt to load specified frame from file (0 = first frame; -1 = last frame)
    int        loadFrame(int f)  { lock(); int r=reader.loadFrame(simul, f); unlock(); return r; }

    /// load next frame in file
    int        loadNextFrame()   { lock(); int r=reader.loadNextFrame(simul); unlock(); return r; }

    
    /// return the Single that is manipulated by the User
    Single const* handle() const;

    /// make a new Single that can be controlled by the user
    Single *   createHandle(Vector const&, real range);
    
    /// switch current handle
    bool       selectClosestHandle(Vector const&, real range);
    
    /// detach current handle
    void       detachHandle();
    
    /// move the current handle
    void       moveHandle(Vector const&);
    
    /// move all handles
    void       moveHandles(Vector const&);
    
    /// delete all handles
    void       deleteHandles();
    
    /// detach current handle from mouse control
    void       releaseHandle() { mHandle = nullptr; }
    
};


#endif
