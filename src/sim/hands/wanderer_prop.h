// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef WANDERER_PROP_H
#define WANDERER_PROP_H

#include "digit_prop.h"

/// Additional Property for Wanderer
/**
 @ingroup Properties
 */
class WandererProp : public DigitProp
{
    friend class Wanderer;
    
public:
    
    /// Diffusion rate on the fiber
    real    diffusion;
    
private:
    
    /// derived variable
    real    diff_rate_dt;
    
    /// derived variable
    real    diff2_rate_dt;
    
    /// derived variable
    real    step_size_kT;
    
    /// derived variable
    real    diff_rate;
    
    /// derived variable
    real    diff_rate_2;
    
    /// taken from the couple, the delta G associated with a stretch of one step size divided by 2kT.
    real    U_step_kT_2;
    
public:
    
    /// constructor
    WandererProp(const std::string& n) : DigitProp(n) { clear(); }
    
    /// destructor
    ~WandererProp() { }
    
    /// return a Hand with this property
    virtual Hand * newHand(HandMonitor* h) const;
    
    /// set default values
    void clear();
    
    /// set from a Glossary
    void read(Glossary&);
    
    /// compute values derived from the parameters
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new WandererProp(*this); }
    
    /// write all values
    void write_values(std::ostream &) const;
    
    /// Complete Hand's property given stiffness value from the Couple/Single
    void completeStiffness(Simul const&, real stiffness);
};

#endif

