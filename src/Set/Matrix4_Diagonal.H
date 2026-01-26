#ifndef SET_MATRIX4_DIAGONAL_H
#define SET_MATRIX4_DIAGONAL_H

#include "AMReX_Config.H"
#include "Util/Util.H"
#include "Base.H"

#include "Set.H"
#include "Set/Matrix3.H"
#include "Set/Matrix4.H"

namespace Set
{
template<>
class Matrix4<AMREX_SPACEDIM, Sym::Diagonal>
{
    std::array<std::array<Set::Scalar,AMREX_SPACEDIM>,AMREX_SPACEDIM> A = {{
            AMREX_D_DECL(
                {AMREX_D_DECL(NAN,NAN,NAN)},
                {AMREX_D_DECL(NAN,NAN,NAN)},
                {AMREX_D_DECL(NAN,NAN,NAN)}
                )}};
        
    
public:
    Set::Matrix getA() const
    {
        Set::Matrix ret;
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                ret(i,j)=A[i][j];
        return ret;
    }

    AMREX_GPU_HOST_DEVICE Matrix4() = default;
    //AMREX_GPU_HOST_DEVICE Matrix4(Set::Matrix a_A): A(a_A) {};
    AMREX_FORCE_INLINE
        Scalar operator () (const int i, const int j, const int k, const int l) const
    {
        if (i == k && j == l) return A[i][j];
        else return 0.0;
    }
    void Randomize()
    {
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                A[i][j] = Util::Random();
    }
    static Matrix4<AMREX_SPACEDIM,Sym::Diagonal> Random()
    {
        Matrix4<AMREX_SPACEDIM,Sym::Diagonal> ret;
        ret.Randomize();
        return ret;
    }
    void Print(std::ostream& os)
    {
        os << "A = \n" << getA();
    }
    static Matrix4<AMREX_SPACEDIM, Sym::Diagonal> Zero()
    {
        Matrix4<AMREX_SPACEDIM, Sym::Diagonal> zero;
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                zero.A[i][j] = 0.0;
        return zero;
    }
    static Matrix4<AMREX_SPACEDIM, Sym::Diagonal> Ones()
    {
        Matrix4<AMREX_SPACEDIM, Sym::Diagonal> ones;
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                ones.A[i][j] = 1.0;
        return ones;
    }
    friend Set::Matrix operator * (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Set::Matrix& b);
    friend Set::Vector operator * (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Set::Matrix3& b);
    friend bool operator == (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& b);
    friend Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator + (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& b);
    friend Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator - (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& b);
    friend Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator * (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Set::Scalar& b);
    friend Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator * (const Set::Scalar& b, const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a);
    friend Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator / (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Set::Scalar& b);


    //AMREX_GPU_HOST_DEVICE void operator =  (Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda =  a.lambda; mu =  a.mu;}
    AMREX_GPU_HOST_DEVICE void operator += (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a) { 
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                A[i][j] += a.A[i][j]; 
    }
    AMREX_GPU_HOST_DEVICE void operator -= (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a) { 
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                A[i][j] -= a.A[i][j]; 
    }
    AMREX_GPU_HOST_DEVICE void operator *= (const Set::Scalar& alpha) {
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                A[i][j] *= alpha; 
    }
    AMREX_GPU_HOST_DEVICE void operator /= (const Set::Scalar& alpha) { 
        for (int i = 0 ; i < AMREX_SPACEDIM; i++)
            for (int j = 0 ; j < AMREX_SPACEDIM; j++)
                A[i][j] /= alpha; 
    }

    Set::Scalar Norm()
    {
        return getA().lpNorm<2>();
    }
    bool contains_nan() const
    {
        if (std::isnan(getA().lpNorm<2>())) return true;
        return false;
    }

};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Set::Matrix operator * (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Set::Matrix& b)
{
    Set::Matrix ret;

    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        for (int j = 0; j < AMREX_SPACEDIM; j++)
        {
            ret(i, j) = a.A[i][j] * b(i, j);
        }
    }

    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Set::Vector operator * (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a, const Set::Matrix3& b)
{
    Set::Vector ret = Set::Vector::Zero();

    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j = 0; j < AMREX_SPACEDIM; j++)
            ret(i) += a.A[i][j] * b(i, j, j);

    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
bool operator == (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a,
    const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& b)
{
    return a.A == b.A;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator + (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a,
    const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& b)
{
    Matrix4<AMREX_SPACEDIM, Sym::Diagonal> ret = a;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j = 0; j < AMREX_SPACEDIM; j++)
            ret.A[i][j] += b.A[i][j];
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator - (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a,
    const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& b)
{
    Matrix4<AMREX_SPACEDIM, Sym::Diagonal> ret = a;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j = 0; j < AMREX_SPACEDIM; j++)
            ret.A[i][j] -= b.A[i][j];
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator * (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a,
    const Set::Scalar& b)
{
    Matrix4<AMREX_SPACEDIM, Sym::Diagonal> ret = a;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j = 0; j < AMREX_SPACEDIM; j++)
            ret.A[i][j] *= b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator * ( const Set::Scalar& b,
                                                    const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a)
{
    return a*b;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Matrix4<AMREX_SPACEDIM, Sym::Diagonal> operator / (const Matrix4<AMREX_SPACEDIM, Sym::Diagonal>& a,
    const Set::Scalar& b)
{
    Matrix4<AMREX_SPACEDIM, Sym::Diagonal> ret = a;
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j = 0; j < AMREX_SPACEDIM; j++)
            ret.A[i][j] /= b;
    return ret;
}


}
#endif
