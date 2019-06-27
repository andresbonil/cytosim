// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef TUBULE_PROP
#define TUBULE_PROP

#include "fiber_prop.h"

class Tubule;


/// additional Property for Tubule
/**
 @ingroup Properties
 */
class TubuleProp : public FiberProp
{
    friend class Tubule;
    
public:
    
    /**
     @defgroup TubulePar Parameters of Tubule
     @ingroup Parameters
     Inherits @ref FiberPar.
     @{
     */
    
    /// see @ref Tubule

    /// Model for dynamic assembly ([0]=PLUS_END, [1]=MINUS_END)
    int     dynamic_model[2];
    
    /// Characteristic force for polymer assembly (default=+inf) ([0]=PLUS_END, [1]=MINUS_END)
    real    growing_force[2];
    
    /// switching rate to the growing state for a fiber shorter than `min_length` (default=0)
    real    rebirth_rate[2];
    
    real    dynamic_speed1[4];
    real    dynamic_speed2[4];
    real    dynamic_trans1[4];
    real    dynamic_trans2[4];
    
    /// @}
    
    /// derived variable
    real    rebirth_prob[2];

private:
    

public:
    
    /// constructor
    TubuleProp(const std::string& n) : FiberProp(n) { clear(); }

    /// destructor
    ~TubuleProp() { }
    
    /// return a Fiber with this property
    Fiber* newFiber() const;
    
    /// set default values
    void clear();
       
    /// set using a Glossary
    void read(Glossary&);
   
    /// check and derive parameter values
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new TubuleProp(*this); }

    /// write
    void write_values(std::ostream &) const;

};

#endif

