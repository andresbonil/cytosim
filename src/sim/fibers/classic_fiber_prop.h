// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef CLASSIC_FIBER_PROP
#define CLASSIC_FIBER_PROP

#include "fiber_prop.h"


/// Enables support for an option to make catastrophe rate dependent on fiber length
#define NEW_LENGTH_DEPENDENT_CATASTROPHE 0

/// Enables support for an option that induces catastrophe if the PLUS_END is outside
#define NEW_CATASTROPHE_OUTSIDE 0

/// Enables support for an option that induces catastrophe following the accumulation of motors in time at the tip of a microtubule following the model of Tischer et al. 2010
#define NEW_CATASTROPHE_TIP_MOTORS 1

/// additional Property for ClassicFiber
/**
 @ingroup Properties
 */
class ClassicFiberProp : public FiberProp
{
    friend class ClassicFiber;
    
public:
    
    /**
     @defgroup ClassicFiberPar Parameters of ClassicFiber
     @ingroup Parameters
     Inherits @ref FiberPar.
     The first column of numbers applies to PLUS_END, and the second to MINUS_END
     @{
     */
    
    /// see @ref ClassicFiber
    
    /// Speed of assembly state
    /**
     Antagonistic force decrease assembly rate exponentially if it is directed against the assembly:

         if ( force < 0 )
             speed = growing_speed * free_polymer * exp( force / growing_force ) + growing_off_speed;
         else
             speed = growing_speed * free_polymer + growing_off_speed;
     
     The parameters are:
     - `growing_speed`, the force-dependent and concentration-dependent assembly rate.
     - `growing_off_speed`, a constant term, normally negative to represent spontaneous disassembly.
     - `growing_force`, the characteristic force
     .
     In this equation, `free_polymer` represents the fraction of free monomers in [0,1].
     Antagonistic force is negative ( force < 0 ) if it is directed against fiber assembly.
     */
    real    growing_speed[2];

    /// Constant term in the growing speed equation
    real    growing_off_speed[2];

    
    /// Characteristic force of assembly state (default=+inf)
    /**
     Antagonistic force decrease assembly rate exponentially.
     */
    real    growing_force[2];
    
    
    /// speed of disassembly state
    /**
     Disassembly occurs always at the specified speed:

         speed = shrinking_speed;

     */
    real    shrinking_speed[2];
    
    
    /// Rate of stochastic switching from assembly to disassembly
    /**
     The catastrophe rate depends on the growth rate of the corresponding tip,
     which is itself reduced by antagonistic force:

         catastrophe_rate_real = catastrophe_rate_stalled / ( 1 + coef * growing_speed_real )

     where `growth_speed_real` is calculated as explained in @ref growing_speed,
     and `coef` is set to match the given `catastrophe_rate` in the absence of force:

         coef = ( catastrophe_rate_stalled/catastrophe_rate - 1.0 ) / growing_speed_unloaded
         growing_speed_unloaded = growing_speed + growing_off_speed;
     
     Note that if `catastrophe_rate_stalled >> catastrophe_rate`, the equation simplies to

         catastrophe_rate_real = catastrophe_rate * growing_speed_unloaded / growing_speed_real

     */
    real    catastrophe_rate[2];

    /// Rate of catastrophe when the growth is stalled
    /**
     If this parameter is not set, the catastrophe rate will not depend on growth speed.
     */
    real    catastrophe_rate_stalled[2];

#if NEW_CATASTROPHE_OUTSIDE
    
    /// catastrophe rate scaling factor applied if the PLUS_END is outside
    /**
     A value < 1 inhibits catastrophe at the edge; A value > 1 accelerates catastrophes
     */
    real    catastrophe_outside;

    /// space used for `catastrophe_outside'
    std::string catastrophe_space;

#endif
    
#if NEW_LENGTH_DEPENDENT_CATASTROPHE
    
    /// Switch to enable the length-dependent catastrophe rate
    /**
     If this is defined, the catastrophe rate will depend on the length of the fiber:

         catastrophe_rate_real = catastrophe_rate * length() / catastrophe_length;

     */
    real    catastrophe_length;
    
#endif
    
#if NEW_CATASTROPHE_TIP_MOTORS
    
    /// Switch to enable the length-dependent catastrophe rate based on the model of Tischer et al. 2010 for the yeast anaphase.
    /**
    If this is defined, the catastrophe rate will depend on the difference of length of the fiber from the moment it started to grow (length()-length_at_rescue). We choose this assumption because the catastrophy rate in the anaphase spindle increases as the microtubule grows, but it seems to depend only on the time the microtubule has grown, and not on the length of the spindle, suggesting that kinesin-8 molecules are recruited either at the midzone, or after the midzone.

        catastrophe_rate_real = catastrophe_rate * (1.-psi*exp(-L/Lm)-(1.-psi)*exp(-L/Lt))
        psi = Lm/(Lm-Lt)
        Lm: (v_m/k_u) -> Run length of the motor: velocity of the motor divided by unbinding rate of the motor
        Lt: v_g/beta -> Growth speed of the microtubule divided by unbinding rate of the motor at the tip.
        catastrophe_rate: Represents the catastrophe rate of a microtubule of infinite length
        
    */
    
    /// Run length of the motor: velocity of the motor divided by unbinding rate of the motor (Lm= v_m/k_u)
    real catastrophe_Lm;
    
    /// Growth speed of the microtubule divided by unbinding rate of the motor at the tip (Lt= v_g/beta)
    real catastrophe_Lt;
    
#endif
    
    /// Rate of stochastic switching from disassembly to assembly
    real    rescue_rate[2];
    
    /// switching rate to the growing state for a fiber shorter than `min_length` (default=0)
    real    rebirth_rate[2];
    
    /// @}
    
private:
    
    real    shrinking_speed_dt[2];
    real    growing_speed_dt[2];
    real    growing_off_speed_dt[2];
    real    catastrophe_rate_dt[2];
    real    catastrophe_rate_stalled_dt[2];
    real    catastrophe_coef[2];
    real    rescue_prob[2], rebirth_prob[2];
    real    catastrophe_psi;
    
#if NEW_CATASTROPHE_OUTSIDE
    /// pointer to actual Space used for `catastrophe_outside`
    Space const* catastrophe_space_ptr;
#endif
    
public:
    
    /// constructor
    ClassicFiberProp(const std::string& n) : FiberProp(n) { clear(); }

    /// destructor
    ~ClassicFiberProp() { }
    
    /// return a Fiber with this property
    Fiber* newFiber() const;
    
    /// set default values
    void clear();
       
    /// set using a Glossary
    void read(Glossary&);
   
    /// check and derive parameter values
    void complete(Simul const&);
    
    /// return a carbon copy of object
    Property* clone() const { return new ClassicFiberProp(*this); }

    /// write
    void write_values(std::ostream&) const;

};

#endif

