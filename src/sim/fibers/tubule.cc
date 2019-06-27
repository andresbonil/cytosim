// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#include "dim.h"
#include "sim.h"
#include "assert_macro.h"
#include "tubule.h"
#include "tubule_prop.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "simul.h"
#include "space.h"


//========================================================================
//  - - - - - - - - - - - - - - CONSTRUCTORS - - - - - - - - - - - - - - -
//========================================================================

void Tubule::resetDynamics()
{
#ifdef OLD_TRACK_CONTACTS
    mtContactTime = -1;
#endif
    
    mtState[ PLUS_END] = STATE_GREEN;
    mtState[MINUS_END] = STATE_WHITE;
    
    mtGrowth[ PLUS_END] = 0;
    mtGrowth[MINUS_END] = 0;
}


Tubule::Tubule(TubuleProp const* p) : Fiber(p), prop(p)
{
    resetDynamics();
}


Tubule::~Tubule()
{
    //Cytosim::log("destroying Tubule %p\n", this);
    prop = nullptr;
}


//------------------------------------------------------------------------------
#pragma mark -


void Tubule::setDynamicStateM(unsigned s)
{
    if ( s!=STATE_WHITE && s!=STATE_GREEN && s!=STATE_RED )
        throw InvalidParameter("Invalid AssemblyState for Tubule MINUS_END");
    
    if ( s != mtState[MINUS_END] )
    {
        mtState[MINUS_END] = s;
    }

}


void Tubule::setDynamicStateP(unsigned s)
{
    if ( s!=STATE_WHITE && s!=STATE_GREEN && s!=STATE_RED )
        throw InvalidParameter("Invalid AssemblyState for Tubule PLUS_END");
    
    if ( s != mtState[PLUS_END] )
    {
        mtState[PLUS_END] = s;
        
#ifdef OLD_TRACK_CONTACTS
        /* instrumentation code, Francois Nedelec, Oct. 2005
         records state transitions, with position and timing */
        
        if ( s==STATE_RED )
        {
            FILE * file = fopen("logs/contacts.txt", "a");
            if ( file )
            {
                if ( !ferror(file) )
                {
                    //report the timing of the catastrophe:
                    real time = simul().time();
                    fprintf(file, "%-4lu %8.3f %8.3f ", identity(), mtContactTime, time);
                    
                    //report position, and direction at PLUS_END
                    posEndP().print(file);
                    dirEndP().println(file);
                }
                fclose(file);
            }
            mtContactTime = -1;
        }
#endif
    }
}

int  Tubule::dynParamIndex(const FiberEnd which) const
{
    assert_true( mtState[which]==STATE_GREEN || mtState[which]==STATE_RED );
    
    if ( which==PLUS_END )
    {
        if ( mtState[which] == STATE_GREEN )
            return 0;
        else
            return 1;
    }
    else
    {
        assert_true( which==MINUS_END );
        if ( mtState[which] == STATE_GREEN )
            return 2;
        else
            return 3;
    }
    
}


real Tubule::givenGrowthRate(const FiberEnd which, const int alt) const
{
    if ( mtState[which] == STATE_WHITE ) 
        return 0;
    
    if ( alt )
        return prop->dynamic_speed2[ dynParamIndex(which) ];
    else
        return prop->dynamic_speed1[ dynParamIndex(which) ];
}

real Tubule::givenTransitionRate(const FiberEnd which, const int alt) const
{
    if ( mtState[which] == STATE_WHITE )
        return 0;
    
    if ( alt )
        return prop->dynamic_trans2[ dynParamIndex(which) ];
    else
        return prop->dynamic_trans1[ dynParamIndex(which) ];
}


/**
 - model = 2:
 transitions depend on the MT-length,
 there is no rescue beyond L=13 microns, and no catastrophy at L=0
 the rates given for rescue/catastrophies are those at L=10 microns.
 Dogterom, Felix, Guet, Leibler, J. Cell Biol. 1996 (mitotic X.Eggs extracts)
 Typical rates for this option should be:  cata=0.03 /sec. resc=0.01 /sec.
 - model = 3:
 Transitions depend on the growth rate of the MT tip, itself reduced by
 antagonistic force. the correspondance is: 1/rate = a + b * growthSpeed
 cf. `Dynamic instability of MTs is regulated by force`
 M.Janson, M. de Dood, M. Dogterom. JCB 2003, Figure 2 C.

 */
real Tubule::transitionRate(FiberEnd which, const int model) const
{
    real growth = mtGrowth[which] / prop->time_step;
    real rate0 = givenTransitionRate(which);
    real rate1 = givenTransitionRate(which, 1);
    
    switch ( model )
    {
        case 0:
            // rate is zero: no transition occurs
            return 0;
        
        case 1:
            // transitions are independent of MT length 
            return rate0;
        
        case 2:
            if ( givenGrowthRate(which) > 0 )
                return rate0 * length() * 0.1;
            else {
                if ( length() < 13.0 )
                    return rate0 * (13.0-length()) * 0.333;
                else
                    return 0;
            }
        break;
        
        
        case 3: {
            // 1 / catastrophe_rate depends linearly on growing speed
            real grow0 = givenGrowthRate(which);
            real grow1 = givenGrowthRate(which, 1);
            real coef = 0;
            if ( rate0 > 0 )
                coef = ( rate1 / rate0 - 1.0 ) / ( grow0 + grow1 );
            return rate1 / ( 1 + coef * growth );
        }
        
        case 7:
            //Transitions are alternate if a Hand is attached
            if ( nbHandsNearEnd(0.1, which) )
                return rate1;
            else
                return rate0;
            
        case 9:
            /** the position of the plus-end in the Space determines its dynamics:
            transitions are trans1[] inside the box and trans2[] outside
            */
            if ( prop->confine_space_ptr->outside(posEnd(which)) )
                return rate1;
            else
                return rate0;
            
            
        default: {
            throw InvalidParameter("invalid value tubule:model=", model);
        }
    }
    return 0;
}


//------------------------------------------------------------------------------
#pragma mark -

void Tubule::stepPlusEnd()
{
    const FiberEnd which = PLUS_END;
    
    if ( mtState[which] != STATE_WHITE )
    {    
        // get force acting on the point at the end:
        real force  = projectedForceEndP();
        real speed0 = givenGrowthRate(which);
        real speed1 = givenGrowthRate(which, 1);
        real speed;
        
        if ( speed0 > 0 )
        {
            // antagonistic force (projectedForceEnd < 0) decrease assembly exponentially
            if ( force < 0  &&  prop->growing_force[0] < INFINITY )
                speed = speed0 * prop->free_polymer * exp(force/prop->growing_force[0]) + speed1;
            else
                speed = speed0 * prop->free_polymer + speed1;
        }
        else
        {
            // disassembly is constant
            speed = speed0 + speed1;
        }
        
        mtGrowth[which] = prop->time_step * speed;
 
        //switch to the other state with the probability given by transitionRate():
        real rate = transitionRate(which, prop->dynamic_model[0]);
        
        if ( RNG.test( rate * prop->time_step ))
            setDynamicStateP((mtState[which]==STATE_GREEN) ? STATE_RED : STATE_GREEN);
    }
}


void Tubule::stepMinusEnd()
{
    const FiberEnd which = MINUS_END;
    
    if ( mtState[which] != STATE_WHITE )
    {    
        // get force acting on the point at the end:
        real force  = projectedForceEndM();
        real speed0 = givenGrowthRate(which);
        real speed1 = givenGrowthRate(which, 1);
        real speed;
        
        if ( speed0 > 0 )
        {
            // antagonistic force (projectedForceEnd < 0) decrease assembly exponentially
            if ( force < 0  &&  prop->growing_force[1] < INFINITY )
                speed = speed0 * prop->free_polymer * exp(force/prop->growing_force[1]) + speed1;
            else
                speed = speed0 * prop->free_polymer + speed1;
        }
        else
        {
            // disassembly is constant
            speed = speed0 + speed1;
        }
        
        mtGrowth[which] = prop->time_step * speed;
        
        //switch to the other state with the probability given by transitionRate():
        real rate = transitionRate(which, prop->dynamic_model[1]);
        
        if ( RNG.test( rate * prop->time_step ))
            setDynamicStateM((mtState[which]==STATE_GREEN) ? STATE_RED : STATE_GREEN);
    }
}


void Tubule::step()
{    
    if ( prop->dynamic_model[0] )
    {
        stepPlusEnd();
        
        if ( length() + mtGrowth[PLUS_END] > prop->min_length )
            growP(mtGrowth[PLUS_END]);
        else
        {
            // do something if the fiber is too short:
            if ( !prop->persistent )
            {
                delete(this);
                // exit to avoid doing anything with a dead object:
                return;
            }
            
            if ( RNG.test(prop->rebirth_prob[0]) )
                setDynamicStateP(STATE_GREEN);
        }
    }
    
    if ( prop->dynamic_model[1] )
    {
        stepMinusEnd();
    
        if ( length() + mtGrowth[MINUS_END] > prop->min_length )
            growM( mtGrowth[MINUS_END] );
        else
        {
            if ( !prop->persistent )
            {
                delete(this);
                // exit to avoid doing anything with a dead object:
                return;
            }
            
            if ( RNG.test(prop->rebirth_prob[1]) )
                setDynamicStateM(STATE_GREEN);
        }
    }

#ifdef OLD_TRACK_CONTACTS
    //intrumentation code for S. pombe, Nedelec Feb 2006
    //initiate contact-time of plus-ends with cell ends:
    if ( isGrowing(PLUS_END)  &&  mtContactTime < 0 )
    {
        //check for contact:
        Space const* spc = prop->confine_space_ptr;
        Vector w = posEndP();
        
        if ( spc->outside(w)  &&  fabs(w.XX) > spc->length(0) )
            mtContactTime = simul().time();
    }
#endif

    Fiber::step();
}
                  
                  
//------------------------------------------------------------------------------
#pragma mark -

void Tubule::write(Outputter& out) const
{
    Fiber::write(out);

    /// write variables describing the dynamic state of the ends:
    writeHeader(out, TAG_DYNAMIC);
    out.writeUInt16(mtState[MINUS_END]);
    out.writeUInt16(mtState[PLUS_END]);
}


void Tubule::read(Inputter & in, Simul& sim, ObjectTag tag)
{
    resetDynamics();

#ifdef BACKWARD_COMPATIBILITY
    if ( tag == TAG_DYNAMIC || ( tag == TAG && in.formatID() < 44 ) )
#else
    if ( tag == TAG_DYNAMIC )
#endif
    {
#ifdef BACKWARD_COMPATIBILITY
        if ( in.formatID() < 42 && in.formatID() > 30 )
        {
            mtState[MINUS_END] = in.readUInt8();
            mtState[PLUS_END]  = in.readUInt8();
        }
        else
#endif
        {
            // current format
            mtState[MINUS_END] = in.readUInt16();
            mtState[PLUS_END]  = in.readUInt16();
        }
    }
#ifdef BACKWARD_COMPATIBILITY
    if ( tag != TAG_DYNAMIC || in.formatID() < 44 )
#else
    else
#endif
        Fiber::read(in, sim, tag);
}

