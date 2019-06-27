// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "dynamic_bundle.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "messages.h"
#include "simul_prop.h"
extern SimulProp & simprop;
#include "random.h"
#include "tubule.h"
#include <cmath>

//------------------------------------------------------------------------------
//this was done for the study of C.elegans spindle dynamics by C.Kozlowski

void TubuleBundle::mtbInit(const int nbMT)
{
    if ( nbMT < 0 )
        throw InvalidParameter("TUBULE_BUNDLE: mtmax[1] should be >= 2");
    if ( nbMT > mtbAllocated )
        throw InvalidParameter("TUBULE_BUNDLE: mtmax[1] is too large: increase mtbAllocated and recompile");
    
    MSG_ONCE("OPTION TUBULE_BUNDLE MTs are bundles\n");
    
    mtbMax = nbMT;
    for ( int ii = 0; ii < mtbAllocated; ++ii ) {
        mtbLength[ii] = SP.mtinitlength;
        mtbState[ii]  = MTB_SEED;
    }
    
    current = 0;
}

int TubuleBundle::mtbNbTubules() const
{
    int result = 0;
    for ( int ii = 0; ii < mtbMax; ++ii )
        if ( mtbState[ii] != MTB_SEED ) ++result;
    return result;
}


void TubuleBundle::mtbStepTubule(MTBState& state, real& length)
{
    //perform basic MT-dynamic simulation for one of the sub-microtubule
    switch( state ) {
        
        case MTB_GROW:
#ifdef ASSEMBLY_STEPS
            length += ASSEMBLY_STEPS * RNG.poisson( prop->dynanmic_speed1[0]*simprop.dt/ASSEMBLY_STEPS );
#else
            length += SP.mtdynspeed1[0] * simprop.dt;
#endif
            if ( SP.mtdyntrans1[0]>0  && RNG.test(SP.mtdyntrans1[0]*simprop.dt) )
                state = MTB_SHRINK;
            break;
            
            
        case MTB_SHRINK:
#ifdef ASSEMBLY_STEPS
            length -= ASSEMBLY_STEPS * RNG.poisson( -SP.mtdynspeed1[1]*simprop.dt/ASSEMBLY_STEPS );
#else
            length += SP.mtdynspeed1[1] * simprop.dt;
#endif
            if ( length < SP.mtminlength ) {
                state = MTB_SEED;
                break;
            }
            if ( SP.mtdyntrans1[1]>0  && RNG.test(SP.mtdyntrans1[1]*simprop.dt) )
                state = MTB_GROW;
            break;
            
            
        case MTB_SEED:
            //rescue with rate SP.mtnucleationrate:
            if ( RNG.test(SP.mtnucleationrate_dt) )
                state = MTB_GROW;
            break;
            
            
        case MTB_PAUSE:
            break;
    }
}


void TubuleBundle::mtbStepDynamic(int& state, real& length)
{
    mtbState[current]  = static_cast<MTBState>(state);
    mtbLength[current] = length;
    
    int next = current;
    //we skip zero, to have exactly mtbMax fibers in the bundle
    for ( int ii = 1; ii < mtbMax; ++ii ) {
        
        if ( ii != current )
            mtbStepTubule(mtbState[ii], mtbLength[ii]);
        
        if ( mtbLength[ii] > length ) {
            next = ii;
            length  = mtbLength[ii];
            state   = mtbState[ii];
        }
    }
    current = next;
}


void TubuleBundle::mtbWrite(OutputWrapper & out, const Name name) const
{
    assert( name > 0 );
    if ( mtbMax <= 0 ) return;
    out.write("\nb");    assert( 'b' == TAG );
    out.writeUInt16(name, false);
    out.writeUInt16(mtbMax);
    for ( int ii = 0; ii < mtbMax ; ++ii ) {
        out.writeUInt8(mtbState[ii]);
        out.writeReal32(mtbLength[ii]);
    }
}


void TubuleBundle::mtbRead(InputWrapper & in)
{
    mtbMax = in.readUInt16();
    for ( int ii = 0; ii < mtbMax ; ++ii ) {
        mtbState[ii] = (MTBState)in.readUInt8();
        mtbLength[ii] = in.readReal32();
    }
}
