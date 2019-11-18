// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TREADMILLING_FIBER_H
#define TREADMILLING_FIBER_H

#include "sim.h"
#include "vector.h"
#include "node_list.h"
#include "fiber.h"

class TreadmillingFiberProp;


/// A Fiber with assembly at both ends 
/**
 This is not documented yet!

 See the @ref TreadmillingFiberPar.

 @todo Document TreadmillingFiber
 @ingroup FiberGroup
 */
class TreadmillingFiber : public Fiber
{   
private:
    
    /// state of PLUS_END
    state_t    mStateP;
    
    /// assembly during last time-step
    real       mGrowthP;
    
    /// state of MINUS_END
    state_t    mStateM;
    
    /// assembly during last time-step
    real       mGrowthM;
    
public:
    
    /// the Property of this object
    TreadmillingFiberProp const* prop;
  
    /// constructor
    TreadmillingFiber(TreadmillingFiberProp const*);

    /// destructor
    virtual ~TreadmillingFiber();
        
    //--------------------------------------------------------------------------
    
    /// return assembly/disassembly state of MINUS_END
    state_t     dynamicStateM() const;
    
    /// change state of MINUS_END
    void        setDynamicStateM(state_t s);
    
    /// the amount of freshly assembled polymer during the last time step
    real        freshAssemblyM() const;

    
    /// return assembly/disassembly state of PLUS_END
    state_t     dynamicStateP() const;

    /// change state of PLUS_END
    void        setDynamicStateP(state_t s);
    
    /// the amount of freshly assembled polymer during the last time step
    real        freshAssemblyP() const;
    
    //--------------------------------------------------------------------------
    
    /// monte-carlo step
    void        step();
    
    //--------------------------------------------------------------------------
    
    /// write to Outputter
    void        write(Outputter&) const;

    /// read from Inputter
    void        read(Inputter&, Simul&, ObjectTag);
    
};


#endif
