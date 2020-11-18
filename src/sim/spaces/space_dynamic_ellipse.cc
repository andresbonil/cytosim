// Cytosim 3.0 - F. Nedelec and Laboratory, Copyright EMBL 2007

#include "dim.h"
#include "space_dynamic_ellipse.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "glossary.h"
#include "smath.h"

/// prefactor for volume computation: Pi in 2D and 4/3 Pi in 3D
constexpr real VPREF = (DIM+1)*M_PI/3.0;

/// power for ellipsoid surface calculation
constexpr real POW = 1.6075;

/// building block for area of an ellipsoid in 3D
inline real surf_block(const real a, const real b)
{
    return pow(a*b,POW);
}

/// building block for area of an ellipsoid in 3D
inline real surf_block(const real a, const real b, const real c)
{
    return pow(a*b,POW) + pow(b*c,POW) + pow(a*c,POW);
}


SpaceDynamicEllipse::SpaceDynamicEllipse(SpaceDynamicEllipseProp const* p)
: SpaceEllipse(p) , prop(p)
{
    if ( DIM == 1 )
        throw InvalidParameter("dynamic_ellipse is not usable in 1D");
    
    pressure = 0;
    mat = MatrixD::identity();
    inv = MatrixD::identity();
    
    reset_forces();
    rad_forces.set(0,0,0);
}

//-------------------------------------------------------------------------------------
//        Set interactions and update forces felt by ellipse.
//-------------------------------------------------------------------------------------

void SpaceDynamicEllipse::setInteractions(Meca&) const
{
    reset_forces();
}


/**
 Call the appropriate interaction from \a meca, to confine \a pe, which is at position \a pos.
 
 The default implementation projects \a pos,
 to calculate the direction of the normal to the edge of the Space,
 and then calls Meca::interPlane, with the approprimate aguments.
 This generates a friction-less potential centered on the edge.

 Also update \a Rforces and \a Torques that will be use to evolve the Space
*/
void SpaceDynamicEllipse::setInteraction(Vector const&pos, Mecapoint const& pe, Meca& meca, real stiff) const
{
    Vector prj;
    prj = project(pos);
    Vector dir = pos - prj;
    real n = dir.normSqr();
    if ( n > 0 )
    {
        // Register the force to the ellipse
        decompose_force(stiff * dir, prj, dir);
        // And to the meca
        meca.addPlaneClamp(pe, prj, dir, stiff/n);
    }
}



//------------------------------------------------------------------------------
//        Computing forces
//------------------------------------------------------------------------------

void SpaceDynamicEllipse::reset_forces() const
{
    Torques = nullTorque;
    Rforces.set(0,0,0);
}

/**
 register forces applied to the space
 */
void SpaceDynamicEllipse::decompose_force(Vector const& forces, Vector const& pos,Vector const& dir) const
{
#if ( 0 )
    // keep only force component in the normal direction:
    Vector nfo = dir * ((forces*dir) / dir.normSqr());
    add_radial_force(nfo, proj);
#else
    add_radial_force(forces, pos);
#endif
    Torques += cross(pos, forces);
}

/**
 Add a point-like force acting on the ellipse
 */
void SpaceDynamicEllipse::add_radial_force(Vector const& forces, Vector const& pos) const
{
    Vector U = director(0);
    Rforces.XX += dot(U, forces) * dot(U, pos) / length_[0];
#if ( DIM >= 2 )
    Vector V = director(1);
    Rforces.YY += dot(V, forces) * dot(V, pos) / length_[1];
#endif
#if ( DIM > 2 )
    Vector W = director(2);
    Rforces.ZZ += dot(W, forces) * dot(W, pos) / length_[2];
#endif
}


// ----------------------------------------------
//  Internal forces
// ----------------------------------------------


/**
 The derivative of pressure energy with respect to each ellipse parameter:
 EP = Volume * Pressure
 dEP/da = dV/da * Pressure
*/
Vector SpaceDynamicEllipse::pressure_forces(const real P) const
{
    const real S = VPREF * P;
#if ( DIM == 1 )
    return Vector(0, 0, 0);
#elif ( DIM == 2 )
    return Vector(S*length_[1], S*length_[0]);
#elif ( DIM > 2 )
    return Vector(S*length_[1]*length_[2], S*length_[2]*length_[0], S*length_[0]*length_[1]);
#endif
}


/*
 Pressure is a Lagrange multiplier associated with volume conservation
 We follow Newtons's method to minimize:
 
     F = Volume(next_time_step) - prop->volume
 
 Hence we iterate:
 
     P = P - F / dF/dP
 
 until the machine precision is exhausted
*/
real SpaceDynamicEllipse::compute_pressure(Vector const& sizes, Vector const& radif) const
{
    real P = pressure;
    real err = INFINITY, last_err;
    
    if ( prop->mobility_dt <= 0 )
        return 0;
    
    size_t cnt = 0;
    do {
        last_err = err;
        
        // the objective is to reach desired volume at the next time-step:
        Vector dim = sizes + ( radif + pressure_forces(P) ) * prop->mobility_dt;
        
        real der = VPREF * VPREF * prop->mobility_dt;
        
        real r0 = dim.XX;
#if ( DIM == 2 )
        // r0 = A0 + P * B0;   B0 = VPREF * r1 * mobility_dt
        // r1 = A1 + P * B1;   B1 = VPREF * r0 * mobility_dt
        // vol = VPREF * ( A0 + P * B0 ) * ( A1 + P * B1 );
        // der = VPREF * ( B0 * r1 + r0 * B1 )
        real r1 = dim.YY;
        err = VPREF * r0 * r1 - prop->volume;
        der *= square(r0) + square(r1);
#elif ( DIM > 2 )
        // r0 = A0 + P * B0;   B0 = VPREF * r1 * r2 * mobility_dt
        // r1 = A1 + P * B1;   B1 = VPREF * r2 * r0 * mobility_dt
        // r2 = A2 + P * B2;   B2 = VPREF * r0 * r1 * mobility_dt
        // vol = VPREF * ( A0 + P * B0 ) * ( A1 + P * B1 ) * ( A2 + P * B2 );
        // der = VPREF * ( B0 * r1 * r2 + r0 * B1 * r2 + r0 * r1 * B2 )
        real r1 = dim.YY;
        real r2 = dim.ZZ;
        err = VPREF * r0 * r1 * r2 - prop->volume;
        der *= square(r0*r1) + square(r0*r2) + square(r1*r2);
#endif

        P -= err / der;
        
        if ( ++cnt > 256 )
        {
            std::clog << "pressure calculation failed at " << err << '\n';
            return 0;
        }

    } while ( fabs(err) < fabs(last_err) );
    //std::clog << "volume error " <<  err << '\n';

    return P;
}


/**
 The derivative of surface energy with respect to each ellipse parameter:
 ES = Surface * tension
 dES/da = dS/da * tension
 */
Vector SpaceDynamicEllipse::tension_forces() const
{
    Vector res;

#if ( DIM == 2 )

    real S = -prop->tension * M_PI;
    real N = sqrt( (3.0*length_[0]+length_[1])*(length_[0]+3.0*length_[1]) );

    res.XX = S * (3.0 - ( 3.0*length_[0] + 5.0*length_[1] ) / N );
    res.YY = S * (3.0 - ( 3.0*length_[1] + 5.0*length_[0] ) / N );
    
#elif ( DIM > 2 )
    
    real S = -prop->tension * surfaceEllipse(length_);
    
    real pXY = surf_block(length_[0], length_[1]);
    real pXZ = surf_block(length_[0], length_[2]);
    real pYZ = surf_block(length_[1], length_[2]);
    real XYZ = surf_block(length_[0], length_[1], length_[2]);

    res.XX = S * ( pXY + pXZ ) / ( length_[0] * XYZ );
    res.YY = S * ( pXY + pYZ ) / ( length_[1] * XYZ );
    res.ZZ = S * ( pXZ + pYZ ) / ( length_[2] * XYZ );

#endif
    return res;
}

//-------------------------------------------------------------------------------------
///    Update ellipse shape
//-------------------------------------------------------------------------------------
    

void SpaceDynamicEllipse::step()
{
    if ( prop->volume > 0 )
    {
        rad_forces = Rforces;
        
        // calculate forces:
        Rforces += tension_forces();
        pressure = compute_pressure(length_, Rforces);
        Rforces += pressure_forces(pressure);

        // implement changes in shape:
        if ( prop->mobility_dt > 0 )
        {
            Vector delta = prop->mobility_dt * Rforces;
            for ( int i=0; i<DIM; ++i )
            {
                assert_true(delta[i] == delta[i]);
                length_[i] += delta[i];
            }
        }
        
        // implement rotation:
        if ( prop->mobility_rot_dt > 0 )
        {
            //std::clog << "% DynamicEllipse torque " << Torques << "\n";
#if ( DIM == 2 )
            real theta = prop->mobility_rot_dt * Torques;
            if ( theta > REAL_EPSILON )
            {
                real c = cos(theta);
                real s = sin(theta);
                mat = Matrix22(c, s, -s, c) * mat;
            }
#elif ( DIM > 2 )
            real n = Torques.norm();
            real theta = prop->mobility_rot_dt * n;
            if ( theta > REAL_EPSILON )
            {
                MatrixD rot = MatrixD::rotationAroundAxis(Torques/n, cos(theta), sin(theta));
                mat = rot * mat;
            }
#endif
            // update rotation:
            inv = mat.transposed();
        }
    }
    
    reset_forces();
}




/// Checking consistency of ellipse sizes
void SpaceDynamicEllipse::resize(Glossary& opt)
{
    SpaceEllipse::resize(opt);
    prop->volume = volume();
    //std::cout << " dynamic_ellipse:volume set to " << prop->volume << std::endl;
    
    if ( prop->volume <= 0 )
    {
        prop->volume = volume();
        //std::cout << " dynamic_ellipse:volume set to " << prop->volume << std::endl;
    }
    
}


//-------------------------------------------------------------
// Utilities
//-------------------------------------------------------------


Vector SpaceDynamicEllipse::director(unsigned ix) const
{
    assert_true(ix < DIM);
    return mat.column(ix);
}


real SpaceDynamicEllipse::surfaceEllipse(Vector const& sizes)
{
#if ( DIM > 2 )
    real a = sizes.XX;
    real b = sizes.YY;
    real c = sizes.ZZ;
    return 4.0*M_PI*pow(surf_block(a,b,c)/3.0, 1.0/POW);
#elif ( DIM == 2 )
    // In 2D, the 'surface' is a line
    real a = sizes.XX;
    real b = sizes.YY;
    real r = square( (a-b)/(a+b) );
    return M_PI*(a+b)*(1.0+3.0*r/(10.0+sqrt(4.0-3.0*r)));
#else
    return 0;
#endif
}


real SpaceDynamicEllipse::volumeEllipse(Vector const& sizes)
{
#if ( DIM > 2 )
    return VPREF*sizes[0]*sizes[1]*sizes[2];
#elif ( DIM == 2 )
    return VPREF*sizes[0]*sizes[1];
#endif
    return 0;
}


void SpaceDynamicEllipse::read(Inputter& in, Simul& sim, ObjectTag tag)
{
    SpaceEllipse::read(in, sim, tag);
    unsigned n = in.readUInt16();
    if ( n != 10 )
        throw InvalidIO("Unexpected data in SpaceDynamicEllipse::read");
    real vol = in.readFloat();
    // adjust volume in Property:
    prop->volume = vol;
    // read 3x3 orientation matrix:
#if ( DIM > 2 )
    for ( unsigned j = 0; j < 3; ++j )
    for ( unsigned i = 0; i < 3; ++i )
        mat(i,j) = in.readFloat();
#else
    real m[9] = { 0 };
    for ( unsigned i = 0; i < 9; ++i )
        m[i] = in.readFloat();
#if ( DIM == 2 )
    mat(0,0) = m[0];
    mat(1,0) = m[1];
    mat(0,1) = m[3];
    mat(1,1) = m[4];
#endif
#endif
}


void SpaceDynamicEllipse::write(Outputter& out) const
{
    SpaceEllipse::write(out);
    out.writeUInt16(10);
    out.writeFloat(prop->volume);
#if ( DIM > 2 )
    // write matrix in column-major format
    for ( unsigned j = 0; j < 3; ++j )
    for ( unsigned i = 0; i < 3; ++i )
        out.writeFloat(mat(i,j));
#else
    real m[9] = { 0 };
#if ( DIM == 2 )
    m[0] = mat(0,0);
    m[1] = mat(1,0);
    m[3] = mat(0,1);
    m[4] = mat(1,1);
#endif
    for ( unsigned i = 0; i < 9; ++i )
        out.writeFloat(m[i]);
#endif
}



//------------------------------------------------------------------------------
//                         OPENGL  DISPLAY
//------------------------------------------------------------------------------

#ifdef DISPLAY

#include "gle.h"

bool SpaceDynamicEllipse::draw() const
{
#if ( 0 )
    // display principal axes:
    glBegin(GL_LINES);
    for ( unsigned n=0; n < DIM; ++n )
    {
        glVertex2f(0,0);
        gle::gleVertex(length(n)*director(n));
    }
    glEnd();
#endif

    real MM[16] = { 0 };
    MM[ 5]=1.0f;
    MM[10]=1.0f;
    MM[15]=1.0f;
    for ( unsigned i=0; i<DIM; ++i )
    for ( unsigned j=0; j<DIM; ++j )
        MM[i+4*j] = mat(i,j);
    glPushMatrix();
    gle::gleMultMatrix(MM);
    SpaceEllipse::draw();
    glPopMatrix();

    return true;
}

#else

bool SpaceDynamicEllipse::draw() const
{
    return false;
}


#endif

