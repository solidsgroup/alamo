#ifndef SET_MATRIX4_ISOTROPIC_H
#define SET_MATRIX4_ISOTROPIC_H

#include "Util/Util.H"
#include "Base.H"

namespace Set
{
template<>
class Matrix4<AMREX_SPACEDIM,Sym::Isotropic>
{
    Set::Scalar lambda=NAN, mu=NAN;
public:
    AMREX_GPU_HOST_DEVICE Matrix4() {};
    AMREX_GPU_HOST_DEVICE Matrix4(Set::Scalar a_lambda, Set::Scalar a_mu) : lambda(a_lambda), mu(a_mu) {};
    /// Note: for the Isotropic Matrix4 this routine works for **retrieval only**!
    /// If you try to assign a value using this with, say
    ///
    ///     isotropicmatrix4(i,j,k,l) = 8.0
    ///
    /// you will get a `lvalue required as left operand of assignment` compile error.
    /// You should probably consider using a lower symmetry operator.
    AMREX_FORCE_INLINE
    Scalar operator () (const int i, const int j, const int k, const int l) const
    {
        Set::Scalar ret = 0.0;
        if (i==k && j==l) ret += mu;
        if (i==l && j==k) ret += mu;
        if (i==j && k==l) ret += lambda;
        return ret;
    }
    void Randomize()
    {
        lambda = Util::Random();
        mu = Util::Random();
    }
    void Print (std::ostream& os )
    {
        os << "lambda = " << lambda << " mu = " << mu;
    }
    Set::Scalar Lambda () const
    {
        return lambda;
    }
    Set::Scalar Mu () const
    {
        return mu;
    }
    static Matrix4<AMREX_SPACEDIM,Sym::Isotropic> Zero()
    {
        Matrix4<AMREX_SPACEDIM,Sym::Isotropic> zero;
        zero.lambda = 0.0;
        zero.mu = 0.0;
        return zero;
    }

    //
    // The isotropic tensor is defined to be
    //
    // C_{ijkl} = mu(d_{il}d_{jk} + d_{ik}d_{jl}) + lambda d_{ij} d_{kl}
    // 
    // The inverse ("compliance") tensor is 
    // 
    // S_{ijkl} = ((1+nu)/2E)(d_{il}d_{jl} + d_{ik}d_{jl}) - (nu/E)d_{ij}d_{kl}
    //
    // Replacing E, nu with Lame constants gives:
    // 
    // S_{ijkl} = (1/(4 mu))(d_{il}d_{jl} + d_{ik}d_{jl}) + (lambda/(2*mu*(3lambda+2mu))) * d_{ij} * d_{kl}
    //
    // For reference: http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm
    //                https://en.wikipedia.org/wiki/Lam%C3%A9_parameters
    //
    Matrix4<AMREX_SPACEDIM,Sym::Isotropic> Inverse()
    {
        Matrix4<AMREX_SPACEDIM,Sym::Isotropic> inv;
        inv.mu = 1./4./mu;
        inv.lambda = lambda / (2*mu*(3*lambda + 2*mu));
        return inv;
    }
    friend Set::Matrix operator * (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Matrix  &b);
    friend Set::Vector operator * (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Matrix3 &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator - (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator * (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Scalar &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator / (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Scalar &b);
    friend Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator * (const Set::Scalar &b, const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a);
    friend Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator + (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a,const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &b);

    //AMREX_GPU_HOST_DEVICE void operator =  (Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda =  a.lambda; mu =  a.mu;}
    AMREX_GPU_HOST_DEVICE void operator += (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda += a.lambda; mu += a.mu;}
    AMREX_GPU_HOST_DEVICE void operator -= (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda -= a.lambda; mu -= a.mu;}
    AMREX_GPU_HOST_DEVICE void operator *= (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda *= a.lambda; mu *= a.mu;}
    AMREX_GPU_HOST_DEVICE void operator /= (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a) {lambda /= a.lambda; mu /= a.mu;}
    AMREX_GPU_HOST_DEVICE void operator *= (const Set::Scalar &alpha) {lambda *= alpha; mu *= alpha;}
    AMREX_GPU_HOST_DEVICE void operator /= (const Set::Scalar &alpha) {lambda /= alpha; mu /= alpha;}
    
    Set::Scalar Norm()
    {
        return std::sqrt(lambda*lambda + mu*mu);
    }

};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Matrix operator * (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Matrix &b)
{
    Set::Matrix ret;
        
    #if   AMREX_SPACEDIM == 2
    ret(0,0) = (a.lambda + 2.*a.mu) * b(0,0) +       a.lambda      *b(1,1);
    ret(1,1) =        a.lambda      * b(0,0) + (a.lambda + 2.*a.mu)*b(1,1);
    ret(0,1) = a.mu*(b(0,1) + b(1,0)); ret(1,0) = ret(0,1);
    
    #elif AMREX_SPACEDIM == 3
    ret(0,0) = (a.lambda + 2.*a.mu) * b(0,0) +       a.lambda      *b(1,1) +       a.lambda      *b(2,2);
    ret(1,1) =        a.lambda      * b(0,0) + (a.lambda + 2.*a.mu)*b(1,1) +       a.lambda      *b(2,2);
    ret(2,2) =        a.lambda      * b(0,0) +       a.lambda      *b(1,1) + (a.lambda + 2.*a.mu)*b(2,2);
    ret(1,2) = a.mu*(b(1,2) + b(2,1)); ret(2,1) = ret(1,2);
    ret(2,0) = a.mu*(b(2,0) + b(0,2)); ret(0,2) = ret(2,0);
    ret(0,1) = a.mu*(b(0,1) + b(1,0)); ret(1,0) = ret(0,1);
    
    #endif         
    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Set::Vector operator * (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Matrix3 &b)
{
    Set::Vector ret = Set::Vector::Zero();
    
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int j=0; j < AMREX_SPACEDIM; j++)
            ret(i) += a.mu*(b(i,j,j) + b(j,i,j)) + a.lambda*b(j,j,i);

    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator * (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Scalar &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Isotropic> ret;// = Set::Vector::Zero();
    ret.mu = a.mu * b;
    ret.lambda = a.lambda * b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator / (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Set::Scalar &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Isotropic> ret;// = Set::Vector::Zero();
    ret.mu = a.mu / b;
    ret.lambda = a.lambda / b;
    return ret;
}
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator * (const Set::Scalar &b, const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a)
{
    return a*b;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator + (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Isotropic> ret;// = Set::Vector::Zero();
    ret.mu = a.mu + b.mu;
    ret.lambda = a.lambda + b.lambda;
    return ret;
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE 
Matrix4<AMREX_SPACEDIM,Sym::Isotropic> operator - (const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &a, 
                                                    const Matrix4<AMREX_SPACEDIM,Sym::Isotropic> &b)
{
    Matrix4<AMREX_SPACEDIM,Sym::Isotropic> ret = a;
    ret.mu -= b.mu;
    ret.lambda -= b.lambda;
    return ret;
}


}
#endif
