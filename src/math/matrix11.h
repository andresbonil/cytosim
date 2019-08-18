// Cytosim was created by Francois Nedelec. Copyright 2007-2017 EMBL.
// F. Nedelec, Strasbourg 08.06.2018

#ifndef MATRIX11
#define MATRIX11

#include "real.h"
#include "vector1.h"
#include <cstdio>
#include <iostream>

/// 1x1 matrix class with 1 'real' elements
/**
 This matrix can use AVX instructions if 'real == double'
 */
class Matrix11
{
public:

    real val;
   
    Matrix11() {}
    
    /// construct Matrix from coordinate
    Matrix11(real a)
    {
        val = a;
    }

    ~Matrix11() {}
    
    /// dimensionality
    static int dimension() { return 1; }
    
    /// leading dimension
    static int stride() { return 1; }

    /// set all elements to zero
    void reset()
    {
        val = 0.;
    }
    
    /// true if element is different from 'zero'
    bool operator != (real zero) const
    {
        return ( val != zero );
    }

    /// copy values from lower triangle to upper triangle
    void copy_lower()
    {
    }

    /// conversion to array of 'real'
    real* data()             { return &val; }
    real* addr(int i, int j) { return &val; }

    /// access functions to element by index
    real& operator[](int i)       { return val; }
    real  operator[](int i) const { return val; }
    
    /// access functions to element by line and column indices
    real& operator()(int i, int j)       { return val; }
    real  operator()(int i, int j) const { return val; }
    
    /// extract column vector at given index
    Vector1 column(const unsigned) const
    {
        return Vector1(val);
    }
    
    /// extract line vector at given index
    Vector1 line(const unsigned) const
    {
        return Vector1(val);
    }

    /// extract diagonal
    Vector1 diagonal() const
    {
        return Vector1(val);
    }

    /// human-friendly output
    void print(FILE * f) const
    {
        fprintf(f, "[ %9.3f ]\n", val);
    }
    
    /// conversion to string
    std::string to_string(int w, int p) const
    {
        std::ostringstream os("[ ");
        os.precision(p);
        os << std::setw(w) << std::fixed << val << " ]";
        return os.str();
    }

    /// true is matrix is symmetric
    bool is_symmetric() const
    {
        return true;
    }
    
    /// scale all elements
    void scale(const real alpha)
    {
        val *= alpha;
    }
    
    /// scaled matrix
    const Matrix11 operator *(const real alpha) const
    {
        return Matrix11(val * alpha);
    }

    /// return sum of two matrices
    const Matrix11 operator +(Matrix11 const& M) const
    {
        return Matrix11( val + M.val );
    }

    /// return difference of two matrices
    const Matrix11 operator -(Matrix11 const& M) const
    {
        return Matrix11( val - M.val );
    }

    /// subtract given matrix
    void operator +=(Matrix11 const& M)
    {
        val += M.val;
    }

    /// add given matrix
    void operator -=(Matrix11 const& M)
    {
        val -= M.val;
    }
    
    /// transpose matrix in place
    void transpose()
    {
    }
    
    /// return transposed matrix
    Matrix11 transposed() const
    {
        return Matrix11(val);
    }
    
    /// maximum of all component's absolute values
    real norm() const
    {
        return fabs(val);
    }

    /// multiplication by a vector: this * V
    const Vector1 vecmul(Vector1 const& V) const
    {
        return Vector1(val * V.XX);
    }
    
    /// multiplication by a vector: this * V
    const Vector1 vecmul(real const* ptr) const
    {
        return Vector1(val * ptr[0]);
    }

    /// matrix-vector multiplication
    friend Vector1 operator * (Matrix11 const& mat, Vector1 const& vec)
    {
        return mat.vecmul(vec);
    }

    /// multiplication by a vector: transpose(M) * V
    const Vector1 trans_vecmul(real const* V) const
    {
        return Vector1(val * V[0]);
    }

    /// multiplication by another matrix: @returns this * B
    const Matrix11 mul(Matrix11 const& B) const
    {
        return Matrix11(val * B.val);
    }
    
    /// matrix-matrix multiplication
    friend Matrix11 operator * (Matrix11 const& mat, Matrix11 const& mut)
    {
        return mat.mul(mut);
    }

    /// multiplication by another matrix: @returns transpose(this) * B
    const Matrix11 trans_mul(Matrix11 const& B) const
    {
        return Matrix11(val * B.val);
    }

    /// add full matrix: this <- this + M
    void add_full(Matrix11 const& M)
    {
        val += M.val;
    }
    
    /// add full matrix: this <- this + alpha * M
    void add_full(const real alpha, Matrix11 const& M)
    {
        val += alpha * M.val;
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + M
    void add_half(Matrix11 const& M)
    {
        val += M.val;
    }
    
    /// add lower triangle of matrix including diagonal: this <- this + alpha * M
    void add_half(const real alpha, Matrix11 const& M)
    {
        val += alpha * M.val;
    }
    
    /// subtract lower triangle of matrix including diagonal: this <- this - M
    void sub_diag(Matrix11 const& M)
    {
        val -= M.val;
    }

    
    /// add all elements of block 'S' to array 'M'
    void addto(real * M, unsigned ldd) const
    {
        M[0] += val;
    }
    
    /// add lower elements of this block to upper triangle of 'M'
    void addto_upper(real * M, unsigned ldd) const
    {
        M[0] += val;
    }
    
    /// add lower elements of this block to both upper and lower triangles of 'M'
    void addto_symm(real * M, unsigned ldd) const
    {
        M[0] += val;
    }
    
    /// add all elements of this block to 'M', with transposition
    void addto_trans(real * M, unsigned ldd) const
    {
        M[0] += val;
    }
    
    /// return `a * Identity`
    static Matrix11 diagonal(real a)
    {
        return Matrix11(a);
    }
    
    /// identity matrix
    static Matrix11 identity()
    {
        return diagonal(1);
    }

    /// set a symmetric matrix: [ dir (x) transpose(dir) ]
    static Matrix11 outerProduct(Vector1 const& dir)
    {
        return Matrix11(dir.XX * dir.XX);
    }

    /// set a symmetric matrix: alpha * [ dir (x) transpose(dir) ]
    static Matrix11 outerProduct(Vector1 const& dir, real alpha)
    {
        return Matrix11(dir.XX * dir.XX * alpha);
    }
    
    /// return outer product: [ dir (x) transpose(vec) ]
    static Matrix11 outerProduct(Vector1 const& dir, Vector1 const& vec)
    {
        return Matrix11(dir.XX * vec.XX);
    }
    
    /// return symmetric matrix block :  dia * I + [ dir (x) dir ] * len
    static Matrix11 offsetOuterProduct(const real dia, Vector1 const& dir, const real len)
    {
        return Matrix11(len * dir.XX * dir.XX + dia);
    }
    
    /// return rotation matrix of angle defined by cosinus and sinus
    static Matrix11 rotation(real c, real s)
    {
        return Matrix11(std::copysign(1, c));
    }

    /// rotation angle
    real rotationAngle() const;

    /// return a rotation that transforms (1,0,0) into (-1,0,0)
    static Matrix11 rotation180();

    /// return a rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix11 rotationToVector(const Vector1& vec);
    
    /// return a random rotation that transforms (1,0,0) into `vec` ( norm(vec) should be > 0 )
    static Matrix11 randomRotationToVector(const Vector1& vec)
    {
        return rotationToVector(vec);
    }
    
    /// a random rotation chosen uniformly
    static Matrix11 randomRotation();

};


/// output operator to std::ostream
inline std::ostream& operator << (std::ostream& os, Matrix11 const& M)
{
    os << "[ " << M.val << " ]";
    return os;
}

#endif

