// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "modulo.h"
#include "dim.h"
#include "exceptions.h"

constexpr int PERIODIC_XYZ = ( 1 << DIM ) - 1;
constexpr int PERIODIC_YZ  = ( 1 << (DIM-1) ) - 1;
constexpr int PERIODIC_X   = 1;


/// enable periodicity in dimension 'd'
void Modulo::enable(int d, real size)
{
    if ( size <= 0 )
        throw InvalidParameter("periodic:length must be > 0");
    mMode |= 1<<d;
    mSize[d] = size;
}


const Vector Modulo::periodicity(int d) const
{
    Vector vec(0,0,0);
    if ( d < DIM && ( mMode & 1<<d ))
        vec[d] = 2 * mSize[d];
    return vec;
}


void Modulo::fold(Vector& vec) const
{
    if ( mMode == PERIODIC_XYZ )
    {
        fold(vec.XX, mSize[0]);
#if ( DIM > 1 )
        fold(vec.YY, mSize[1]);
#endif
#if ( DIM > 2 )
        fold(vec.ZZ, mSize[2]);
#endif
    }
    else if ( mMode == PERIODIC_YZ )
    {
#if ( DIM > 1 )
        fold(vec.XX, mSize[0]);
#endif
#if ( DIM > 2 )
        fold(vec.YY, mSize[1]);
#endif
    }
    else if ( mMode == PERIODIC_X )
    {
        fold(vec.XX, mSize[0]);
    }
    else
    {
        if ( mMode & 1 ) fold(vec.XX, mSize[0]);
#if ( DIM > 1 )
        if ( mMode & 2 ) fold(vec.YY, mSize[1]);
#endif
#if ( DIM > 2 )
        if ( mMode & 4 ) fold(vec.ZZ, mSize[2]);
#endif
    }
}


//this makes modulo around the center 'ref'
void Modulo::fold(Vector & pos, Vector const& ref) const
{
    pos -= ref;
    fold(pos);
    pos += ref;
}


//calculate the offset from the canonical image to actual 'pos'
const Vector Modulo::offset(Vector const& pos) const
{
    Vector img = pos;
    fold(img);
    return pos - img;
}


//calculate the canonical image of 'pos' and return the associated shift
void Modulo::foldOffset(Vector& pos, Vector& off) const
{
    Vector vec = pos;
    fold(pos);
    off = vec - pos;
}
