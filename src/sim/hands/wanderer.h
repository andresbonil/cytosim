// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef WANDERER_H
#define WANDERER_H

#include "digit.h"

class WandererProp;

/**
 The Wanderer class was created to model Ase1/PRC1 diffusible crosslinkers
 by Manuel Lera-Ramirez in 2018--2019
 */
class Wanderer : public Digit
{
private:
    
    /// disabled default constructor
    Wanderer();
    
    /// Gillespie countdown timer for stepping
    real   nextStep;
    
    /// Calculate propensities (k+ = k0 * p_plus, k- = k0 * p_minus)
    void   calcPropensities(Vector const& force, real& p_plus, real& p_minus);

public:
    
    /// Property
    WandererProp const* prop;
    
    /// constructor
    Wanderer(WandererProp const* p, HandMonitor* h);
    
    /// destructor
    ~Wanderer() {}
    
    /// attach and update variables
    void   attach(FiberSite const& site);
    
    /// simulate when `this` is attached but not under load
    void   stepUnloaded();
    
    /// simulate when `this` is attached and under load
    void   stepLoaded(Vector const& force, real force_norm);
    
};

#endif

