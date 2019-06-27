// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.

#include "matrix33.h"
#include "vector2.h"
#include "random.h"

Vector3 Matrix33::rotationAxis() const
{
    return Vector3(val[5]-val[7], val[6]-val[2], val[1]-val[3]);
}

real Matrix33::rotationAngle() const
{
    real trace = val[0] + val[4] + val[8];
    return acos(0.5*(1-trace));
}


void Matrix33::getEulerAngles(real& a, real& b, real& c) const
{
    real cb = sqrt(val[0] * val[0] + val[1] * val[1]);
    
    b = atan2( -val[2], cb );
    
    if ( cb != 0 ) {
        a = atan2( val[1], val[0] );
        c = atan2( val[5], val[8] );
    }
    else {
        a = 0;
        c = atan2( -val[3], val[4] );
    }
}


/// returns a rotation of angle PI around axis Z
Matrix33 Matrix33::rotation180()
{
    return Matrix33(-1, 0, 0, 0, -1, 0, 0, 0, 1);
}


Matrix33 Matrix33::rotationAroundX(const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    return Matrix33(1, 0, 0, 0, c, s, 0, -s, c);
}

Matrix33 Matrix33::rotationAroundY(const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    return Matrix33(c, 0, -s, 0, 1, 0, s, 0, c);
}

Matrix33 Matrix33::rotationAroundZ(const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    return Matrix33(c, s, 0, -s, c, 0, 0, 0, 1);
}


Matrix33 Matrix33::rotationAroundPrincipalAxis(unsigned ii, const real angle)
{
    real c = cos(angle);
    real s = sin(angle);
    
    ii %= 3;
    int jj = (ii+1)%3;
    int kk = (ii+2)%3;
    
    Matrix33 res;
    res.reset();
    res.val[ii+BLD*ii] = 1;
    res.val[jj+BLD*jj] = c;
    res.val[kk+BLD*kk] = c;
    res.val[jj+BLD*kk] = -s;
    res.val[kk+BLD*jj] = s;
    return res;
}


Matrix33 Matrix33::rotationFromAngles(const real a[3])
{
    real ca = cos(a[0]), sa = sin(a[0]);
    real cb = cos(a[1]), sb = sin(a[1]);
    real cc = cos(a[2]), sc = sin(a[2]);
    
    Matrix33 res;

    res.val[0+BLD*0] =  ca*cb;
    res.val[1+BLD*0] =  sa*cb;
    res.val[2+BLD*0] = -sb;
    
    res.val[0+BLD*1] =  ca*sb*sc - sa*cc;
    res.val[1+BLD*1] =  sa*sb*sc + ca*cc;
    res.val[2+BLD*1] =  cb*sc;
    
    res.val[0+BLD*2] =  ca*sb*cc + sa*sc;
    res.val[1+BLD*2] =  sa*sb*cc - ca*sc;
    res.val[2+BLD*2] =  cb*cc;
    
    return res;
}


Matrix33 Matrix33::rotationAroundAxisEuler(const real a[3])
{
    real ca = cos(a[0]), sa = sin(a[0]), ca1 = 1 - ca;
    real cb = cos(a[1]), sb = sin(a[1]);
    real cc = cos(a[2]), sc = sin(a[2]);
    
    real sacc        = sa * cc,           sasc        = sa * sc;
    real saccsb      = sacc * sb,         sacccb      = sacc * cb;
    real ccccca1     = cc * cc * ca1,     ccscca1     = cc * sc * ca1;
    real sbccscca1   = sb * ccscca1,      cbccscca1   = cb * ccscca1;
    real cbcbccccca1 = cb * cb * ccccca1, cbsbccccca1 = cb * sb * ccccca1;
    
    Matrix33 res;
    res.val[0+BLD*0] =  cbcbccccca1 + ca;
    res.val[0+BLD*1] =  cbsbccccca1 - sasc;
    res.val[0+BLD*2] =  cbccscca1   + saccsb;
    
    res.val[1+BLD*0] =  cbsbccccca1 + sasc;
    res.val[1+BLD*1] =  ca - cbcbccccca1 + ccccca1;
    res.val[1+BLD*2] =  sbccscca1   - sacccb;
    
    res.val[2+BLD*0] =  cbccscca1 - saccsb;
    res.val[2+BLD*1] =  sbccscca1 + sacccb;
    res.val[2+BLD*2] =  1 - ccccca1;
    
    return res;
}


Matrix33 Matrix33::randomRotation()
{
    //James Arvo, Fast random rotation matrices. in Graphics Gems 3.
    real u2 = M_PI * RNG.sreal();
    real u3 = RNG.preal();
    Vector3 V( cos(u2)*sqrt(u3), sin(u2)*sqrt(u3), sqrt(1-u3) );
    return householder(V) * rotationAroundZ(M_PI*RNG.sreal());
}


Matrix33 Matrix33::randomRotation(real angle)
{
    return rotationAroundAxis(Vector3::randU(), cos(angle), sin(angle));
}


Matrix33 Matrix33::rotationToVector(const Vector3& vec)
{
    Matrix33 res;
    Vector3 X, Y, Z = normalize(vec);
    Z.orthonormal(X, Y);
    res.setColumns(Z, X, Y);
    return res;
}


Matrix33 Matrix33::randomRotationToVector(const Vector3& vec)
{
    Matrix33 res;
    Vector3 X, Y, Z = normalize(vec);
    Z.orthonormal(X, Y);
#if ( 0 )
    real a = M_PI * RNG.sreal();
    real c = cos(a), s = sin(a);
    res.setColumns(Z, X*c+Y*s, Y*c-X*s);
#else
    Vector2 cs = Vector2::randU();
    res.setColumns(Z, X*cs.XX+Y*cs.YY, Y*cs.XX-X*cs.YY);
#endif
    return res;
}

