// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "fiber_segment.h"
#include "space.h"
#include "fiber.h"
#include "simul.h"
#include "modulo.h"

extern Modulo const* modulo;

//------------------------------------------------------------------------------
//---------------- DISTANCE FROM A POINT TO A SECTION OF FIBER -----------------
//------------------------------------------------------------------------------

/**
 W is projected on the line that supports this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the distance between W and its projection
 .
 
 It is assumed here that len() returns the distance between the two points of the FiberSegment
 Attention: `dis` is not set if ( abs < 0 ) or ( abs > len() )
 */
real FiberSegment::projectPoint0(Vector const& w, real& dis) const
{
    assert_true( fib_ );
    
    Vector dx = dir();
    Vector aw = w - pos1();
    
    if ( modulo )
        modulo->fold(aw);
    
    // project with the scalar product:
    real abs = dot(aw, dx);
    
    // calculate distance to projection that fall inside the segment
    if ( 0 <= abs  &&  abs <= len() )
    {
#if ( DIM == 1 )
        dis = 0;
#else
        dis = aw.normSqr() - abs * abs;
#endif
    }
    return abs;
}


/**
 W is projected on the line that supports this FiberSegment
 The function calculates:
 - abs <- the signed distance from pos1() to the projection of W
 - dis <- the distance between W and its projection
 .
 
 It is assumed here that len() returns the distance between the two points of the FiberSegment
 Attention: `dis` may NOT be set if ( abs < 0 ) or ( abs > len() )
 */
real FiberSegment::projectPoint(Vector const& w, real& dis) const
{
    assert_true( fib_ );
    
    Vector dx = dir();
    Vector aw = w - pos1();
    
    if ( modulo )
        modulo->fold(aw);
    
    // project with the scalar product:
    real abs = dot(aw, dx);
    
    // test boundaries of filament:
    if ( abs < 0 )
    {
        if ( isFirst() )
            dis = distanceSqr(w, pos1());
    }
    else if ( abs > len() )
    {
        if ( isLast() )
            dis = distanceSqr(w, pos2());
    }
    else
    {
#if ( DIM == 1 )
        dis = 0;
#else
        dis = aw.normSqr() - abs * abs;
#endif
    }
    return abs;
}


/**
 This may be faster than projectPoint(), but it does not work with periodic boundaries
 */
real FiberSegment::projectPointF(const real w[], real& dis) const
{
    assert_true( fib_ );
    assert_true( !modulo );
    
    const real * p = fib_->addrPoint(pti_);
    
    real dX = p[DIM  ] - p[0];
    real aX = w[0]     - p[0];
#if ( DIM > 1 )
    real dY = p[DIM+1] - p[1];
    real aY = w[1]     - p[1];
#endif
#if ( DIM > 2 )
    real dZ = p[DIM+2] - p[2];
    real aZ = w[2]     - p[2];
#endif
    
    const real ls = len();
    
    // project with the scalar product:
#if ( DIM == 1 )
    real abs = ( dX * aX ) / ls;
#elif ( DIM == 2 )
    real abs = ( dX * aX + dY * aY ) / ls;
#elif ( DIM == 3 )
    real abs = ( dX * aX + dY * aY + dZ * aZ ) / ls;
#endif
    
    // test boundaries of segment:
    if ( abs < 0 )
    {
        if ( isFirst() )
            dis = distanceSqr(w, pos1());
    }
    else if ( abs > ls )
    {
        if ( isLast() )
            dis = distanceSqr(w, pos2());
    }
    else
    {
#if   ( DIM == 1 )
        dis = 0;
#elif ( DIM == 2 )
        dis = aX * aX + aY * aY - abs * abs;
#elif ( DIM == 3 )
        dis = aX * aX + aY * aY + aZ * aZ - abs * abs;
#endif
        
#if ( 0 )
        // verify that the results are identical to projectPoint()
        real d = dis;
        real a = projectPoint(Vector(w), d);
        assert_small(a-abs);
        assert_small(d-dis);
#endif
    }
    return abs;
}


//------------------------------------------------------------------------------
//---------------- DISTANCE TO ANOTHER SECTION OF A FIBER ----------------------
//------------------------------------------------------------------------------

/**
 Evaluate the minimal distance between two segments.
 
 This finds the positions P1, P2 connecting the two supporting lines in the shortest
 way possible. P1 belongs to *this segment, and P2 to that segment. These points 
 are characterized by real abscissa (abs1, abs2) for which values {0, segment_length}
 match the edges of the segments.
 
 The function then evaluates if the points P1, P2 are inside or outside their 
 respective segments.
 
 @return `0` if `P1` or `P2` are outside their respective segment.
 If P1 is outside, then `abs2` is not evaluated.
 
 @return `1` if P1 and P2 are both inside their segments.
 In this case, `dis` is set to be the square of the distance (P1, P2)

 Hence `dis`, `abs1` and `abs2` are only set valid if the return value is `1`.
 If the segments are parallel, the mid-point of the overlapping section is returned.
 */

int FiberSegment::shortestDistance(FiberSegment const& seg, real& abs1, real& abs2, real& dis) const
{
    Vector d1  = diff();
    Vector d2  = seg.diff();
    Vector d12 = seg.pos1() - pos1();
    
    real len1 = len();
    real len2 = seg.len();
    
    if ( modulo )
        modulo->fold(d12);
    
    real beta = dot(d1, d2) / ( len1 * len2 );
    real scal = 1.0 - beta * beta;
    
    if ( scal > REAL_EPSILON )
    {
        // This is the general case of non-parallel lines:
        
        real d1d12 = dot(d1, d12) / ( scal * len1 );
        real d2d12 = dot(d2, d12) / ( scal * len2 );
        
        //abs1 = (( d1 / len1 - beta * d2 / len2 ) * d12 ) / scal;
        abs1 = d1d12 - beta * d2d12;
        
        if ( abs1 < 0  || len1 <= abs1 )
            return 0;
        
        //abs2 = (( beta * d1 / len1 - d2 / len2 ) * d12 ) / scal;
        abs2 = beta * d1d12 - d2d12;
        
        if ( abs2 < 0  || len2 <= abs2 )
            return 0;
        
#if ( DIM > 2 )
        //dis = ( d12 + (abs2/len2) * d2 - (abs1/len1) * d1 ).normSqr();
        //dis = ( d12 - d1 * (abs1/len1) ).normSqr() - abs2 * abs2;
        dis = ( d12 + d2 * (abs2/len2) ).normSqr() - abs1 * abs1;
#else
        dis = 0;
#endif
    }
    else
    {
        /*
         This deals with the case where the two segments are almost parallel:
         beta = +/- 1
         p1 = projection of that.pos1() on this segment
         p2 = projection of that.pos2()
         */
        real p1 = dot(d12, d1) / len1;
        real p2 = p1 + beta * len2;
        
        if ( p1 < 0  &&  p2 < 0 )
            return false;
        
        if ( p1 > len1  &&  p2 > len1 )
            return false;
        
        const real p = p1;
        dis = d12.normSqr() - p1 * p1;
        
        if ( p1 < 0 ) p1 = 0; else if ( p1 > len1 ) p1 = len1;
        if ( p2 < 0 ) p2 = 0; else if ( p2 > len1 ) p2 = len1;
        
        // take the midpoint:
        abs1 = 0.5 * ( p1 + p2 );
        //abs2 = abs1 * beta - ( d12 * d2 ) / len2;
        abs2 = ( abs1 - p ) * beta;
    }
    
    return 1;
}


void FiberSegment::print(std::ostream& os) const
{
    if ( fiber() )
        //os << "(" << fiber()->reference() << " seg " << point() << ":" << point()+1 << ")";
        os << "(f" << fiber()->identity() << " " << point() << ":" << point()+1 << ")";
    else
        os << "(null)";
}


std::ostream& operator << (std::ostream& os, FiberSegment const& obj)
{
    obj.print(os);
    return os;
}
