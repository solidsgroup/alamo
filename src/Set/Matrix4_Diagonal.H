#ifndef SET_MATRIX4_DIAGONAL_H
#define SET_MATRIX4_DIAGONAL_H

#include "Util/Util.H"
#include "Base.H"

namespace Set
{
template<>
class Matrix4<AMREX_SPACEDIM,Sym::Diagonal>
{
    Set::Matrix A = Set::Matrix::Ones();
public:
    AMREX_GPU_HOST_DEVICE Matrix4() {};
    AMREX_GPU_HOST_DEVICE Matrix4(Set::Matrix a_A) : A(a_A) {};
    AMREX_FORCE_INLINE
    Scalar operator () (const int i, const int j, const int k, const int l) const
    {
        if (i==k && j==l) return A(i,j);
        else return 0.0;
    }
    void Randomize()
    {
        A = Set::Matrix::Random();
    }
    void Print (std::ostream& os )
    {
        os << "A = \n" << A;
    }
    static Matrix4<AMREX_SPACEDIM,Sym::Diagonal> Zero()
    {
        Matrix4<AMREX_SPACEDIM,Sym::Diagonal> zero;
        zero.A = Set::Matrix::Zero();
        return zero;
    }
    friend Set::Matrix operator * (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Set::Matrix  &b);
    friend Set::Vector operator * (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Set::Matrix3 &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator + (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator - (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator * (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Set::Scalar &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator / (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Set::Scalar &b);


    //AMREX_GPU_HOST_DEVICE void operator =  (Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda =  a.lambda; mu =  a.mu;}
    AMREX_GPU_HOST_DEVICE void operator += (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a) {A += a.A;}
    AMREX_GPU_HOST_DEVICE void operator -= (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a) {A -= a.A;}
    //AMREX_GPU_HOST_DEVICE void operator *= (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a) {A *= a.A;}
    //AMREX_GPU_HOST_DEVICE void operator /= (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a) {A /= a.A;}
    AMREX_GPU_HOST_DEVICE void operator *= (const Set::Scalar &alpha) {A *= alpha; }
    AMREX_GPU_HOST_DEVICE void operator /= (const Set::Scalar &alpha) {A /= alpha; }

    Set::Scalar Norm()
    {
        return A.lpNorm<2>();
    }

};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Matrix operator * (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Set::Matrix &b)
{
    Set::Matrix ret;
        
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        for (int j = 0; j < AMREX_SPACEDIM; j++)
        {
            ret(i,j) = a.A(i,j) * b(i,j);
        }
    }

    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Vector operator * (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, const Set::Matrix3 &b)
{
    Set::Vector ret = Set::Vector::Zero();
    
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j = 0; j < AMREX_SPACEDIM; j++)
            ret(i) += a.A(i,j) * b(i,j,j);

    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator + (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, 
                                                const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Diagonal> ret = a;
    ret.A += b.A;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator - (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, 
                                                const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Diagonal> ret = a;
    ret.A -= b.A;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator * (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, 
                                                const Set::Scalar &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Diagonal> ret = a;
    ret.A *= b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Diagonal> operator / (const Matrix4<AMREX_SPACEDIM,Sym::Diagonal> &a, 
                                                const Set::Scalar &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Diagonal> ret = a;
    ret.A /= b;
    return ret;
}


}
#endif
