// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
#ifndef RASTERIZER_H
#define RASTERIZER_H

#include "real.h"
#include "vector1.h"
#include "vector2.h"
#include "vector3.h"

/// 2D and 3D rasterizer
/**
 The different functions in the rasterizer call a given method func() for every point 
 of INTEGER coordinates inside a certain volume.
 
 The volume can be specified in two ways:
 - as a polygon described by a list of points, using paintPolygon?D(polygon) 
 - as a cylinder specified by two points P,Q and a scalar 'radius'.
 In 3D, the functions do not rasterize a cylinder, but rather a rectangular volume
 that contains all the points located at distance radius or less from [PQ].

 Important note:   The points defining the polygon do not need to be integers.

 F.Nedelec, EMBL 2002-2017, nedelec@embl.de

*/
namespace Rasterizer 
{
    
    /// a point in 3D, with a bitfield for connectivity
    struct Vertex3
    {
        /// coordinates of the point
        real XX, YY, ZZ;
        
        /// bit-field used to describe the connectivity between the points.
        /**
         Two points A and B are connected if ( A.UU & B.UU ) using the bit-wise AND.
         With a long integer, this limits the number of edges to 64.
         A bigger integer could be used if needed.
        */
        unsigned long UU;
        
        void set(Vector3 const& vec, unsigned long u)
        {
            XX = vec.XX;
            YY = vec.YY;
            ZZ = vec.ZZ;
            UU = u;
        }
        
        void print(std::ostream& os) const
        {
            os << "( " << std::setw(9) << XX;
            os << "  " << std::setw(9) << YY;
            os << "  " << std::setw(9) << ZZ << " )";
        }
    };
    
    /// helper function for qsort() of Vertex3
    int compareVertex3(const void * a, const void * b);

    /// a point in 2D, with slopes in Z
    struct Vertex2
    {
        /// coordinates of the point
        real XX, YY;
        
        /// derivatives of the line with respect to Z
        real dX, dY;
        
        void set(real x, real y, real dx, real dy)
        {
            XX = x;
            YY = y;
            dX = dx;
            dY = dy;
        }
        
        void move()
        {
            XX += dX;
            YY += dY;
        }
        
        void print(std::ostream& os) const
        {
            os << "( " << std::setw(9) << XX;
            os << "  " << std::setw(9) << YY << " )";
        }
    };

    /// accessory data swap function
    template<typename VEC> void swap(VEC& A, VEC& B) { VEC T = A; A = B; B = T; }
    
    /// Compute the Convex Hull of a set of points
    /**
     on entry, pts[] contains 'n_pts' points with coordinates members XX and YY.
     on exit, pts[] is the anti-clockwise polygon hull, starting with the point of lowest Y.
     
     @returns The number of points in the convex hull (at most n_pts).
     */
    template<typename VEC>
    unsigned convexHull2D(unsigned n_pts, VEC pts[])
    {
        //---------- find bottom and top points:
        unsigned inx = 0, top = 0;
        real y_bot = pts[0].YY;
        real y_top = pts[0].YY;
        
        for ( unsigned n = 1; n < n_pts; ++n )
        {
            if ( pts[n].YY < y_bot || ( pts[n].YY == y_bot  &&  pts[n].XX > pts[inx].XX ) )
            {
                inx = n;
                y_bot = pts[n].YY;
            }
            if ( pts[n].YY > y_top || ( pts[n].YY == y_top  &&  pts[n].XX < pts[top].XX ) )
            {
                top = n;
                y_top = pts[n].YY;
            }
        }
        
        if ( inx == top )  //all points are equal ?
            return 1;
        
        // put the bottom point at index zero:
        if ( inx )
        {
            swap(pts[0], pts[inx]);
            if ( top == 0 )
                top = inx;
        }
        // reset
        inx = 0;
        
        // wrap upward on the right side of the hull
        unsigned nxt;
        while ( 1 )
        {
            real pX = pts[inx].XX;
            real pY = pts[inx].YY;
            ++inx;
            
            nxt = top;
            real dx = pts[top].XX - pX;
            real dy = pts[top].YY - pY;
            
            for ( unsigned n = inx; n < n_pts; ++n )
            {
                real dxt = pts[n].XX - pX;
                real dyt = pts[n].YY - pY;
                // keep if slope is lower:
                if ( dxt * dy > dyt * dx )
                {
                    nxt = n;
                    dx = dxt;
                    dy = dyt;
                }
            }
            
            // if we reached bottom point, sweep is complete
            if ( nxt == 0 )
                break;
            
            swap(pts[inx], pts[nxt]);
            
            // if we reached topmost point, change reference for downward sweep:
            if ( nxt == top )
                top = 0;
            else if ( inx == top )
                top = nxt;  // this compensates the swap
        }
        
        return inx;
    }
    
    //------------------------------------ 1D --------------------------------------
#pragma mark - 1D

    /// Rasterizer function in 1D
    void paintFatLine1D(void (*paint)(int, int, int, int, void*),
                        void * arg,            ///< last argument to paint()
                        const Vector1& P,      ///< first end of the segment
                        const Vector1& Q,      ///< second end of the segment
                        real radius,           ///< width by which the segment [PQ] is inflated
                        const Vector1& offset, ///< phase of the grid
                        const Vector1& delta   ///< period for the grid
    );
    
    
    //------------------------------------ 2D --------------------------------------
#pragma mark - 2D
    
    /// Paint a polygon in 2D
    /**
     paintPolygon2D() calls paintPoint(x,y,zz) for every point (x,y) of
     integral coordinates, which are inside the polygon given in xy[]. 
     The polygon should be convex, and ordered anti-clockwise.
     */  
    void paintPolygon2D(void (*paint)(int, int, int, int, void*),
                        void * arg,           ///< last argument to paint
                        unsigned n_pts,       ///< number of points
                        const Vector2[],      ///< the 2D points ( x0 y0 x1 y1 ...)
                        int zz                ///< third coordinate, passed as argument to paint()
                        );
    
    
    /// Paint the inside of a rectangle with edges parallel to the segment PQ
    void paintFatLine2D(void (*paint)(int, int, int, int, void*),
                        void * arg,           ///< last argument to paint
                        const Vector2& P,     ///< first end of the segment [dim=3]
                        const Vector2& Q,     ///< second end of the segment [dim=3]
                        real  length,         ///< length of segment PQ
                        real  radius          ///< width by which PQ is inflated
                        );
    
    
    /// Paint the inside of a rectangle with edges parallel to the segment PQ
    void paintFatLine2D(void (*paint)(int, int, int, int, void*),
                        void * arg,           ///< last argument to paint
                        const Vector2& P,     ///< first end of the segment
                        const Vector2& Q,     ///< second end of the segment
                        real  length,         ///< length of segment PQ
                        const real radius,    ///< width by which the line PQ is extended, to make a round cylinder
                        const Vector2& offset,  ///< phase of the grid
                        const Vector2& delta    ///< period for the grid
                        );
    
    /// Paint a 2D rectangular volume with edges parallel to the main axes
    /**
     The painted volume is square and aligned with the principal axes (X, Y, Z)
     It contains all the points at a distance 'radius' or less from the segment [p,q].
     This is the fastest rasterizer, but the volume can be much greater than that of the cylinder.
     However, the volume is nearly optimal if PQ is aligned with one of the main axis, 
     and paintBox3D is then the best choice.
     */
    void paintBox2D(void (*paint)(int, int, int, int, void*),
                    void * arg,            ///< last argument to paint
                    const Vector2& P,      ///< first end of the segment
                    const Vector2& Q,      ///< second end of the segment
                    real radius,           ///< radius of cylinder contained in the volume
                    const Vector2& offset, ///< phase of the grid
                    const Vector2& delta   ///< period for the grid
                    );
    
    
    //------------------------------------ 3D --------------------------------------
#pragma mark - 3D
    
    /// Polygon rasterizer function in 2D, for Vertex2
    void paintPolygon2D(void (*paint)(int, int, int, int, void*),
                        void * arg,           ///< last argument to paint
                        unsigned nbpts,       ///< number of points
                        const Vertex2[],      ///< points containing additional data
                        int zz = 0            ///< third coordinate, passed as argument to paint()
                        );
    
         
    /// Paint a 3D polygon for which the edges of the convex hull are known
    /**
     The polygon is the convex hull of the 'nbpts' vertices given in pts[].
     Each Vertex contains coordinates and information on the connectivity to other points.
     The connections between Vertices are the edge of the 3D polygon.
     */
    void paintPolygon3D(void (*paint)(int, int, int, int, void*),
                        void * arg,           ///< last argument to paint
                        unsigned n_pts,       ///< number of points
                        Vertex3 pts[]         ///< coordinates + connectivity
                        );
    
    
    /// Paint a 3D cylinder with square section, aligned with the segment [P,Q]
    /**
     A volume is painted around the segment [p,q], containing the cylinder of
     all the points located at a distance 'radius' or less from [p,q].
     The volume is a right cylinder with a square section.
     */
    void paintFatLine3D(void (*paint)(int, int, int, int, void*),
                        void * arg,            ///< last argument to paint
                        const Vector3& P,      ///< first end of the segment
                        const Vector3& Q,      ///< second end of the segment
                        real  length,          ///< length of segment PQ
                        real  radius,          ///< radius of cylinder contained in the volume
                        const Vector3& offset, ///< phase of the grid
                        const Vector3& delta   ///< period for the grid
                        );
    
    
    /// Paint a 3D cylinder with hexagonal section, aligned with the segment [P,Q]
    /**
     A volume is painted around points [p,q], which contains the cylinder of
     all the points at a distance 'radius' or less from the segment [p,q].
     The volume is a right cylinder with hexagonal section.
     This is a tighter approximation of the cylinder than the square cylinder of paintFatLine3D.
     */
    void paintHexLine3D(void (*paint)(int, int, int, int, void*),
                        void * arg,            ///< last argument to paint
                        const Vector3& P,      ///< first end of the segment
                        const Vector3& Q,      ///< second end of the segment
                        real  length,          ///< length of segment PQ
                        real radius,           ///< radius of cylinder contained in the volume
                        const Vector3& offset, ///< phase of the grid
                        const Vector3& delta   ///< period for the grid
                        );

    
    /// Paint a 3D rectangular volume with edges parallel to the main axes
    /**
     The painted volume is square and its edges are parallel to the principal axes (X, Y, Z)
     It contains all the points at a distance 'radius' or less from the segment [p,q].
     This is the fastest rasterizer, but the volume can be much greater than that of the cylinder,
     in particular in the case where PQ >> radius, and PQ is oriented along a diagonal.
     However, the volume is nearly optimal if PQ is almost aligned with one of the main axis, 
     and paintBox3D is then a better choice than the rasterizers that paint a cylinder,
     because it is much faster.
     */
    void paintBox3D(void (*paint)(int, int, int, int, void*),
                    void * arg,           ///< last argument to paint
                    const Vector3& P,     ///< first end of the segment
                    const Vector3& Q,     ///< second end of the segment
                    real radius,          ///< radius of cylinder contained in the volume
                    const Vector3& offset, ///< phase of the grid
                    const Vector3& delta   ///< period for the grid
                   );

}

#endif

