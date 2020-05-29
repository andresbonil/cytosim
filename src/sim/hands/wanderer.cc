// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "digit.h"
#include "wanderer.h"
#include "wanderer_prop.h"
#include "glossary.h"
#include "lattice.h"
#include "simul.h"
extern Random RNG;


Wanderer::Wanderer(WandererProp const* p, HandMonitor* h)
: Digit(p,h), nextStep(0), prop(p)
{
    
}


void Wanderer::attach(FiberSite const& site)
{
    Digit::attach(site);
    nextStep = RNG.exponential();
}


void Wanderer::stepUnloaded()
{
    assert_true( attached() );
    
    // If both positions are occupied, skip
    // We check first whether it is outside, otherwise vacant returns an error.
    int check_m = outsideMP(site()-1) || vacant(site()-1);
    int check_p = outsideMP(site()+1) || vacant(site()+1);

    if ( check_m || check_p )
    {
        nextStep -= prop->diff_rate_dt * ( check_m + check_p );
        
        if ( nextStep <= 0 )
        {
            assert_true( attached() );
            assert_true( check_m + check_p > 0 );
            // -1 or +1, depending on the propensities:
            int step = check_m * check_p * RNG.flipsign() + check_p - check_m;
            
            lati_t s = site() + step;
            
            if ( outsideMP(s) )
            {
                if ( RNG.test_not(prop->hold_growing_end) )
                {
                    detach();
                    return;
                }
            }
            // Here we dont have to check the vacancy of the site anymore, we have checked it already
            else
                hop(s);
            
            nextStep = RNG.exponential();
        }
        
    }
    testDetachment();
}


void Wanderer::calcPropensities(Vector const& force, real& p_plus, real& p_min)
{
    
    // From Hannabuss et al. 2019:
    // step_size_kT = a/kT
    // prop->U_step_kT_2 = (stiffness*a^2)/(2kT)
    // dG0 = f*a/kT
    
    real dG0 = dot(force,dirFiber()) * prop->step_size_kT;
    
    if ( p_plus > 0 )
    {
        real dG = dG0 - prop->U_step_kT_2;
        // To avoid zero / zero
        if ( dG==0 )
            p_plus = 1;
        else
            p_plus = -dG / std::expm1(-dG);
    }
    
    if ( p_min > 0 )
    {
        real dG = -dG0 - prop->U_step_kT_2;
        // To avoid zero / zero
        if ( dG==0 )
            p_min = 1;
        else
            p_min = -dG / std::expm1(-dG);
    }
}


void Wanderer::stepLoaded(Vector const& force, real force_norm)
{
    assert_true( attached() );
    
    
    // If both positions are occupied, skip
    // We check first whether it is outside, otherwise vacant returns an error.
    int check_m = outsideMP(site()-1) || vacant(site()-1);
    int check_p = outsideMP(site()+1) || vacant(site()+1);
    
    if ( check_m || check_p )
    {
        // calculate displacement, dependent on the load along the desired direction of displacement
        real p_plus = 1*check_p;
        real p_min  = 1*check_m;
        
        calcPropensities(force, p_plus, p_min);
        
        nextStep -= ( p_plus + p_min ) * prop->diff_rate_dt;
        
        // I changed this from while to if because otherwise it would move several times when the force is high,
        // otherwise it should update the force every time?
        // This was leading to a lot of steps per dt I guess.
        
        //        if ( (p_plus + p_min)*prop->diff_rate_dt > RNG.pfloat())
        if ( nextStep < 0 )
        {
            assert_true( attached() );
            // -1 or +1, depending on the propensities:
            int step = RNG.flipsign( p_plus/(p_plus+p_min) );
            
            lati_t s = site() + step;
            
            if ( outsideMP(s) )
            {
                if ( RNG.test_not(prop->hold_growing_end) )
                {
                    detach();
                    return;
                }
            }
            // Here we dont have to check the vacancy of the site anymore, we have checked it already
            else
                hop(s);
            nextStep = RNG.exponential();
        }
    }
    
    if ( prop->unbinding_force_inv > 0 )
        testKramersDetachment(force.norm());
    else
        testDetachment();
}

