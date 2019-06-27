// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef TUBULE_BUNDLE
#define TUBULE_BUNDLE

#include "real.h"
#include "object.h"
#include "iowrapper.h"

///should match AssemblyState in tubule.h
enum MTBState
{
    MTB_GROW    = 0,   ///<  Assembly
    MTB_SHRINK  = 1,   ///<  Disassembly
    MTB_PAUSE   = 2,   ///<  Fixed length
    MTB_SEED    = 3    ///<  Fixed length, waiting to nucleate
};

/// a bundle of  dynamic Microtubules
/** This was implemented for the study of C.elegans done by Cleo Kozlowksi,
and published in Kozlowski et al. Cell, May 2007 */
class TubuleBundle
{
private:
    
    ///maximum possible number of sub-MTs in the bundle
    static const int mtbAllocated = 32;
    
    ///actual number of sub-MTs in the bundle
    int            mtbMax;
    
    ///the index of the longest Tubule
    int            current;
    
    ///length of each sub-MT
    real           mtbLength[mtbAllocated];
    
    ///state of each sub-MT
    MTBState       mtbState[mtbAllocated];
    
    ///perform a step to simulate the dynamic of a sub-MT
    void           mtbStepTubule(MTBState&, real&);
        
public:

    ///we need a virtual destructor, since Tubule is derived from TubuleBundle
    virtual       ~TubuleBundle() { }
    
    ///initialization function
    void           mtbInit(int nbMT);
    
    ///perform all the sub-MT steps, and swap to the longest
    void           mtbStepDynamic(int&, real&);
    
  
    ///number of sub-MTs in the bundle
    int            mtbNbSeeds() const { return mtbMax; }
    
    ///length of the ii-th MT in the bundle
    real           mtbTubeLength(const int ii) const { return mtbLength[ii]; }

    ///number of sub-MTs which are not MTB_SEED
    int            mtbNbTubules()  const;
     
    ///state of the ii-th sub-MT plus-end
    MTBState       mtbStateEnd(const int ii) const { return mtbState[ii]; }
    
    //------------------------------ read/write --------------------------------

    /// a unique character identifying the class
    static const Tag TAG = 'b';
    
    /// return unique character identifying the class
    Tag            tag() const { return TAG; }
    
    ///write the sub-MT state to file IO
    void           mtbWrite(OutputWrapper&, const Name) const;
    
    ///write the sub-MT state to file IO
    void           mtbRead(InputWrapper&);

};

#endif
