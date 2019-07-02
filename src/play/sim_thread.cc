// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstdio>
#include <time.h>
#include "sim_thread.h"
#include "exceptions.h"
#include "picket.h"
#include "glapp.h"

//------------------------------------------------------------------------------

/**
 This uses a Parser that cannot write to disc.
 The function callback is called when Parser::hold() is reached.
 */
SimThread::SimThread(void (*callback)(void))
: Parser(simul, 1, 1, 1, 1, 0), hold_callback(callback)
{
    mFlag   = 0;
    mState  = 0;
    mHold   = 0;
    mPeriod = 1;
    pthread_mutex_init(&mMutex, nullptr);
    pthread_cond_init(&mCondition, nullptr);
}

/**
 Possible issue:
 When quitting the application, this destructor might be called after the
 destructor of Simul(), in which case it will access non-existent data,
 most likely causing a crash().
 */
SimThread::~SimThread()
{
    //std::cerr << "~SimThread()\n";
    stop();
    pthread_cond_destroy(&mCondition);
    pthread_mutex_destroy(&mMutex);
    pthread_detach(mChild);
}

//------------------------------------------------------------------------------
#pragma mark - Process control


void SimThread::debug(const char* msg) const
{
    if ( isChild() )
        fprintf(stdout, "- - -   %s\n", msg);
    else
        fprintf(stdout, "+ + +   %s\n", msg);
}


void SimThread::hold()
{
    assert_true( isChild() );

    if ( mFlag )
        terminate();
    
    if ( ++mHold >= mPeriod )
    {
        mHold = 0;
        //debug("holding");
        hold_callback();
        if ( mFlag )
            terminate();
        wait();
    }
}


void SimThread::terminate()
{
    assert_true( isChild() );
    // gently terminate:
    mState = 0;
    unlock();
    //debug("terminate");
    pthread_exit(nullptr);
}


void SimThread::run()
{
    assert_true( isChild() );
    mState = 1;
    
    try {
        if ( Parser::readConfig() )
            std::cerr << "You must specify a config file\n";
    }
    catch( Exception & e ) {
        simul.relax();
        std::cerr << std::endl << "Error: " << e.what() << std::endl;
        //flashText("Error: simulation died");
    }
    // the thread terminates normally:
    hold_callback();
    mState = 0;
    unlock();
}


void* run_helper(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';
    SimThread * lt = static_cast<SimThread*>(arg);
    lt->run();
    return nullptr;
}


void SimThread::extend_run()
{
    assert_true( isChild() );
    mState = 1;
    
    try {
        Parser::execute_run(100000);
    }
    catch( Exception & e ) {
        std::cerr << std::endl << "Error: " << e.what() << '\n';
        simul.relax();
        //flashText("Error: %s", e.what());
    }
    hold_callback();
    mState = 0;
    unlock();
}


void* extend_helper(void * arg)
{
    //std::clog << "slave  " << pthread_self() << '\n';
    SimThread * lt = static_cast<SimThread*>(arg);
    lt->extend_run();
    return nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Process initiation / termination

/**
 This attempts to start the live simulation by
 calling run() in the slave thread
 */
void SimThread::start()
{
    assert_false( isChild() );
    if ( 0 == mState && 0 == trylock() )
    {
        mFlag = 0;
        //std::clog << "master " << pthread_self() << '\n';
        if ( pthread_create(&mChild, nullptr, run_helper, this) )
            throw Exception("failed to create thread");
    }
}


/// call extend_code() in the slave thread
int SimThread::extend()
{
    assert_false( isChild() );
    if ( 0 == mState && 0 == trylock() )
    {
        mFlag = 0;
        //std::clog << "master " << pthread_self() << '\n';
        if ( pthread_create(&mChild, nullptr, extend_helper, this) )
            throw Exception("failed to create thread");
        return 0;
    }
    return 1;
}


void SimThread::step()
{
    assert_false( isChild() );
    if ( mState )
        signal();
}


/**
 ask the slave thread to exit at the next spontaneous halt
*/ 
void SimThread::stop()
{
    assert_false( isChild() );
    if ( mState )
    {
        // request clean termination:
        mFlag = 1;
        signal();
        //debug("waiting on join()");
        // wait for termination:
        pthread_join(mChild, nullptr);
        mState = 0;
    }
}

/**
 kill the slave thread immediately
 */
void SimThread::cancel()
{
    assert_false( isChild() );
    if ( mState )
    {
        mFlag = 2;
        // force termination:
        if ( 0 == pthread_cancel(mChild) )
        {
            // wait for termination:
            pthread_join(mChild, nullptr);
            mState = 0;
            unlock();
        }
    }
}


void SimThread::restart()
{
    assert_false( isChild() );
    stop();
    clear();
    start();
}

//------------------------------------------------------------------------------
#pragma mark - Mouse-controlled Single


SingleProp * SimThread::getHandleProperty() const
{
    Property * p = simul.properties.find("single", "user_single");
    return static_cast<SingleProp*>(p);
}


SingleProp * SimThread::makeHandleProperty(real range)
{
    // Create a Hand that attaches fast and never detach:
    HandProp * hap = new HandProp("user_hand");
    hap->binding_range   = range;
    hap->binding_rate    = 10000;
    hap->unbinding_rate  = 0;
    hap->unbinding_force = INFINITY;
    hap->complete(simul);
    simul.properties.deposit(hap);

    SingleProp * sip = new SingleProp("user_single");
    sip->hand = "user_hand";
    sip->stiffness = 256;
    sip->complete(simul);
    simul.properties.deposit(sip);
    
    return sip;
}


Single * SimThread::createHandle(Vector const& pos, real range)
{
    SingleProp * sip = getHandleProperty();
    if ( !sip )
        sip = makeHandleProperty(range);
    Single * res = new Picket(sip, pos);
    simul.singles.add(res);
    mHandle = res;
    return res;
}


ObjectList SimThread::allHandles(SingleProp const* sip) const
{
    return simul.singles.collect(match_property, sip);
}


bool SimThread::selectClosestHandle(Vector const& pos, real range)
{
    SingleProp * sip = getHandleProperty();
    
    if ( sip )
    {
        real dsm = 0;
        Single * res = nullptr;
        for ( Object * i : allHandles(sip) )
        {
            Single * s = static_cast<Single*>(i);
            real d = ( s->posFoot() - pos ).normSqr();
            if ( !res || d < dsm )
            {
                res = s;
                dsm = d;
            }
        }
        if ( res && dsm < range )
        {
            mHandle = res;
            return 1;
        }
    }
    return 0;
}


Single const* SimThread::handle() const
{
    SingleProp * sip = getHandleProperty();
    if ( sip && mHandle )
    {
        for ( Object * i : allHandles(sip) )
            if ( i == mHandle )
                return mHandle;
    }
    mHandle = nullptr;
    return nullptr;
}


void SimThread::detachHandle()
{
    if ( mHandle )
    {
        if ( mHandle->attached() )
            mHandle->detach();
    }
}

void SimThread::moveHandle(Vector const& pos)
{
    if ( mHandle )
    {
        mHandle->setPosition(pos);
    }
}


void SimThread::moveHandles(Vector const& vec)
{
    SingleProp * sip = getHandleProperty();
    if ( sip )
        ObjectSet::translateObjects(allHandles(sip), vec);
}


void SimThread::deleteHandles()
{
    lock();
    SingleProp * sip = getHandleProperty();
    if ( sip )
        simul.erase(allHandles(sip));
    mHandle = nullptr;
    unlock();
}

void SimThread::clear()
{
    assert_false( isChild() );
    simul.erase();
    mHandle = nullptr;
}

//------------------------------------------------------------------------------
#pragma mark - Parameter modifications

#if ( 0 )

#include <fcntl.h>

/// set file to not block on read() even if data is not available:
void set_nonblocking(int fd)
{
    fcntl(fd, F_SETFL, O_NONBLOCK);
}

#endif


/// check if file has data for input
inline int has_input(int fd)
{
    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(fd, &fds);
    struct timeval tv = {0, 10};   // seconds, microseconds
    return select(1, &fds, nullptr, nullptr, &tv);
}


/**
 Read standard input and executes incoming commands.
 This should be executed by a process who already owns the lock on the data
 */
size_t SimThread::readInput(size_t max_nb_lines)
{
    const size_t LINESIZE = 2048;
    clearerr(stdin);
    
    if ( has_input(STDIN_FILENO) > 0 )
    {
        size_t cnt = 0;
        // some input is available, process line-by-line:
        Parser parser(simul, 1, 1, 1, 1, 1);
        char str[LINESIZE];
        
        // read one line from standard input (including terminating \n):
        while ( fgets(str, LINESIZE, stdin) )
        {
            //write(STDOUT_FILENO, ">>>> ", 5); write(STDOUT_FILENO, str, strlen(str));
            std::stringstream iss(str);
            try {
                parser.evaluate(iss);
                glApp::flashText0(str);
            }
            catch ( Exception & e ) {
                std::cerr << "Error reading stdin: " << e.what() << '\n';
            }
            if ( ++cnt >= max_nb_lines )
                break;
            // check if more input is available:
            if ( has_input(STDIN_FILENO) < 1 )
            {
                //printf("processed %i lines from standard input\n", cnt);
                break;
            }
        }
        return cnt;
    }
    return 0;
}

/**
 Read config file from the start, allowing parameters to be changed, while 
 simulation objects remain as they are. This will pause a running simulation 
 is running live, read the config file, and allow it to proceed.
 */
void SimThread::reloadParameters(std::string const& file)
{
    lock();
    // set a parser that can only change properties:
    if ( Parser(simul, 0, 1, 0, 0, 0).readConfig(file) )
        std::cerr << "Error : File not found";
    //std::cerr << "reloaded " << simul.prop->config_file << std::endl;
    unlock();
}


/**
 This will execute the code specified in `iss`, with full rights to modify Simul.
 
 A simulation running live will be paused; the code is executed in another Parser,
 and the simulation is allowed to proceed.
 */
int SimThread::execute(std::string const& code, std::string const& msg)
{
    lock();
    try {
        std::istringstream iss(code);
        Parser(simul, 1, 1, 1, 1, 1).evaluate(iss);
    }
    catch( Exception & e ) {
        std::cerr << "Error : " << e.what();
        std::cerr << msg;
    }
    unlock();
    return 0;
}

/**
 Save current state in two files
 */
void SimThread::exportObjects(bool binary)
{
    lock();
    try {
        char str[64] = { '\0' };
        
        snprintf(str, sizeof(str), "properties%04i.cmo", reader.currFrame());
        simul.writeProperties(str, true);
        
        snprintf(str, sizeof(str), "objects%04i.cmo", reader.currFrame());
        simul.writeObjects(str, false, binary);
    }
    catch( Exception & e ) {
        std::cerr << "Error in Simul::exportObjects(): " << e.what();
    }
    unlock();
}


void SimThread::writeProperties(std::ostream& os, bool prune)
{
    lock();
    try {
        simul.writeProperties(os, prune);
    }
    catch( Exception & e ) {
        std::cerr << "Error in Simul::writeProperties(): " << e.what();
    }
    unlock();
}

