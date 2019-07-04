// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#ifndef COUPLE_H
#define COUPLE_H

#include "object.h"
#include "hand_monitor.h"
#include "couple_prop.h"
#include "hand.h"

class Meca;


/// A set of two Hand linked by an elastic element
/**
 A Couple contains two pointers to Hand:
 - cHand1
 - cHand2
 .
 There are 4 possible states for a Couple:
 - state FF (0): cHand1 and cHand2 are free,
 - state AF (1): cHand1 is bound, cHand2 is free,
 - state FA (2): cHand1 is free, cHand2 is bound,
 - state AA (3): both hands are attached
 .
 The method state() return the state of the Couple in [0-3].

 Generally the Couple behaves according to its state:
 - FF     : the Couple is diffusing and both Hands are trying to bind fibers,
 - AF, FA : the localization is given by the attachement point on the fiber,
 - AA     : the Couple is acting as a Hookean spring between the two fibers.
 .
 
 The default Couple has:
 - a zero resting length (it uses Meca:interLink())
 - no specificity
 .

 @ingroup CoupleGroup
 */
class Couple : public Object, public HandMonitor
{
public:

    /// associated properties
    CoupleProp const* prop;
    
protected:
    
    /// position and position in previous step of complex
    Vector   cPos;
    
    /// first Hand
    Hand    * cHand1;

    /// second Hand
    Hand    * cHand2;
    
    /// specialization of HandMonitor
    bool      allowAttachment(FiberSite const&);
    /// specialization of HandMonitor
    void      afterAttachment(Hand const*);
    /// specialization of HandMonitor
    void      beforeDetachment(Hand const*);
    /// specialization of HandMonitor
    ObjectID  nucleatorID() const { return Object::identity(); }
    /// Simul container
    Simul*    simul_ptr() const { return &Object::simul(); }
    /// specialization of HandMonitor
    Hand *    otherHand(Hand const*) const;
    /// specialization of HandMonitor
    Vector    otherPosition(Hand const*) const;
    /// specialization of HandMonitor
    Vector    otherDirection(Hand const*) const;
    /// specialization of HandMonitor
    real      interactionLength() const { return prop->length; }
    /// stiffness of the interaction
    real      interactionStiffness() const { return prop->stiffness; }

    /// update position to account for diffusion in one time step
    void      diffuse() { cPos.addRand(prop->diffusion_dt); }
    
public:
    
    /// create following the specifications in the CoupleProp
    Couple(CoupleProp const*, Vector const & w = Vector(0,0,0));

    /// destructor
    virtual ~Couple();

    /// copy operator
    Couple&  operator=(Couple const&);
    
    //--------------------------------------------------------------------------
    
    /// change the property and update the two Hands
    void           setProperty(CoupleProp *);
    
    /// add interactions to the Meca
    virtual void   setInteractions(Meca &) const;
    
    //--------------------------------------------------------------------------
    
    /// the position of the complex, calculated from cPos, cHand1 and cHand2
    virtual Vector position() const;
   
    /// Couple can be displaced only if it is not attached
    virtual int    mobile()               const { return !cHand1->attached() && !cHand2->attached(); }
    
    /// translate object's position by the given vector
    virtual void   translate(Vector const& x)   { cPos += x; }
    
    /// move object to specified position
    void           setPosition(Vector const& x) { cPos = x; }

    /// modulo the current position vector in the space
    virtual void   foldPosition(Modulo const*);
    
    /// set the position randomly inside prop->confine_space
    void           randomizePosition();
    
    //--------------------------------------------------------------------------
    
    /// activity flag
    virtual bool   active()               const { return true; }
    
    /// true if both Hands are attached
    bool           linking()              const { return cHand1->attached() && cHand2->attached(); }

    /// the state of the Couple in { 0 ... 3 } representing { FF, FA, FA, AA }
    int            state()                const { return cHand1->attached() + 2 * cHand2->attached(); }
    
    ///stiffness of the link ( = prop->stiffness )
    real           stiffness()            const;
    
    /// return one of the Hand that is attached, or zero if both are detached
    Hand *         attachedHand()         const;
    
    /// force between hands, essentially: stiffness * ( cHand2->posHand() - cHand1->posHand() )
    virtual Vector force()                const;
     
    /// cosine of the angle between the two Fibers attached by the hands
    real           cosAngle()             const { return dot(cHand1->dirFiber(), cHand2->dirFiber()); }
   
    /// position on the side of fiber1 used for sideInteractions
    virtual Vector posSide()              const { return cHand1->pos(); }
    
    /// the position of the complex if it is unattached
    Vector         posFree()              const { return cPos; }
   
    //--------------------------------------------------------------------------

    /// simulation step for a free Couple: diffusion
    virtual void   stepFF(const FiberGrid&);
    
    /// simulation step for a Couple attached by Hand1
    virtual void   stepAF(const FiberGrid&);
    
    /// simulation step for a Couple attached by Hand2
    virtual void   stepFA(const FiberGrid&);
    
    /// simulation step for a doubly-attached Couple
    virtual void   stepAA();

    //--------------------------------------------------------------------------

    /// pointer to Hand1
    Hand*    hand1()                            { return cHand1; }
    
    /// true if Hand1 is attached
    bool     attached1()                  const { return cHand1->attached(); }
    
    /// Fiber to which Hand1 is attached, or zero if not attached
    Fiber*   fiber1()                     const { return cHand1->fiber(); }
    
    /// attachment position of Hand1 along fiber (only valid if Hand1 is attached)
    real     abscissa1()                  const { return cHand1->abscissa(); }
    
    /// position of Hand1 when attached (only valid if Hand1 is attached)
    Vector   posHand1()                   const { return cHand1->pos(); }
    
    /// direction of Fiber at attachment point of Hand1 (only valid if Hand1 is attached)
    Vector   dirFiber1()                  const { return cHand1->dirFiber(); }
 
    /// attach Hand1 at the given FiberSite
    void     attach1(FiberSite s)               { if ( cHand1->attachmentAllowed(s) ) cHand1->attach(s); }
    
    /// attach Hand1 at the given end
    void     attachEnd1(Fiber* f, FiberEnd end) { cHand1->attachEnd(f, end); }
    
    /// move Hand1 to given end
    void     moveToEnd1(FiberEnd end)           { cHand1->moveToEnd(end); }

    //--------------------------------------------------------------------------

    /// pointer to Hand2
    Hand*    hand2()                            { return cHand2; }
    
    /// true if Hand2 is attached
    bool     attached2()                  const { return cHand2->attached(); }
    
    /// Fiber to which Hand2 is attached, or zero if not attached
    Fiber*   fiber2()                     const { return cHand2->fiber(); }
    
    /// attachment position of Hand2 along fiber (only valid if Hand2 is attached)
    real     abscissa2()                  const { return cHand2->abscissa(); }
    
    /// position of Hand2 when attached (only valid if Hand2 is attached)
    Vector   posHand2()                   const { return cHand2->pos(); }
    
    /// direction of Fiber at attachment point of Hand2 (only valid if Hand2 is attached)
    Vector   dirFiber2()                  const { return cHand2->dirFiber(); }
    
    /// attach Hand2 at the given FiberSite
    void     attach2(FiberSite s)               { if ( cHand2->attachmentAllowed(s) ) cHand2->attach(s); }
    
    /// attach Hand2 at the given end
    void     attachEnd2(Fiber *f, FiberEnd end) { cHand2->attachEnd(f, end); }
    
    /// move Hand2 to given end
    void     moveToEnd2(FiberEnd end)           { cHand2->moveToEnd(end); }

    //--------------------------------------------------------------------------

    /// a static_cast<> of Node::next()
    inline Couple * next()        const { return static_cast<Couple*>(nNext); }
    
    /// a static_cast<> of Node::prev()
    inline Couple * prev()        const { return static_cast<Couple*>(nPrev); }
    
    //------------------------------ read/write --------------------------------

    /// a unique character identifying the class
    static const ObjectTag TAG = 'c';
    
    /// return unique character identifying the class
    ObjectTag       tag() const { return TAG; }
    
    /// return associated Property
    Property const* property() const { return prop; }
    
    /// write to file
    void            write(Outputter&) const;
    
    /// read from file
    void            read(Inputter&, Simul&, ObjectTag);
    
    /// return PointDisp of Hand1
    PointDisp const* disp1() const { return cHand1->prop->disp; }
    
    /// return PointDisp of Hand2
    PointDisp const* disp2() const { return cHand2->prop->disp; }
    
    /// return PointDisp of Hand1
    PointDisp const* disp12() const;
    
    /// return PointDisp of Hand2
    PointDisp const* disp21() const;

};


#endif

