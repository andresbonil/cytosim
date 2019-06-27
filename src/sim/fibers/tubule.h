// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TUBULE_H
#define TUBULE_H

#include "sim.h"
#include "vector.h"
#include "node_list.h"
#include "common.h"
#include "object.h"
#include "fiber.h"


class Single;
class TubuleProp;


// this option records some statistics on contact of plus-end with the cortex
//#define OLD_TRACK_CONTACTS

/// Adds dynamic instability (growth/shrinkage of the ends) to a Fiber
/**
 This code is outdated, and we discourage you to use it.
 @ingroup FiberGroup
 */
class Tubule : public Fiber
{
private:
    
    /// state of FiberEnd
    unsigned   mtState[3];
    
    /// amount growth by the ends (in length-units) during the time-step
    real       mtGrowth[3];
        
    /// resets the values of member variables
    void       resetDynamics();
    
#ifdef OLD_TRACK_CONTACTS
    ///instrumentation code for S. pombe
    real       mtContactTime;
#endif

public:
    
    /// the Property of this object
    TubuleProp const* prop;
        
    /// returns the index which should be used to get parameters
    int dynParamIndex(const FiberEnd which) const;
    
    //--------------------------------------------------------------------------
  
    /// constructor
    Tubule(TubuleProp const*);

    /// destructor
    virtual ~Tubule();
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of the end `which`
    unsigned    dynamicStateM() const { return mtState[MINUS_END]; }

    /// change state of MINUS_END
    void        setDynamicStateM(unsigned s);
    
    /// the amount of freshly assembled polymer during the last time step
    real        freshAssemblyM() const { return mtGrowth[MINUS_END]; }
    
    
    /// return assembly/disassembly state of the end `which`
    unsigned    dynamicStateP() const { return mtState[PLUS_END]; }
    
    /// change state of PLUS_END
    void        setDynamicStateP(unsigned s);
    
    /// the amount of freshly assembled polymer during the last time step
    real        freshAssemblyP() const { return mtGrowth[PLUS_END]; }
    
    //--------------------------------------------------------------------------
    
    /// the growth rate set for the specified end
    real        givenGrowthRate(FiberEnd which, int alt=0) const;
    
    /// the transition rate set for the specified end
    real        givenTransitionRate(FiberEnd which, int alt=0) const;
    
    /// calculate the transition rate of a MT end, function of the parameter, growth and position
    real        transitionRate(FiberEnd which, int model) const;
    
    //--------------------------------------------------------------------------
    
    /// simulate dynamic instability at MINUS_END
    void        stepMinusEnd();
    
    /// simulate dynamic instability at PLUS_END
    void        stepPlusEnd();
    
    /// monte-carlo step
    void        step();
    
    //--------------------------------------------------------------------------
    
    /// write to Outputter
    void        write(Outputter&) const;

    /// read from Inputter
    void        read(Inputter&, Simul&, ObjectTag);
    
};


#endif
