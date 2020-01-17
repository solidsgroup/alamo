#ifndef SET_MATRIX4_ISOTROPIC_H
#define SET_MATRIX4_ISOTROPIC_H

#include "Util/Util.H"
#include "Base.H"

namespace Set
{
template<>
class Matrix4<3,Sym::Isotropic>
{
    Set::Scalar lambda=NAN, mu=NAN;
public:
    AMREX_GPU_HOST_DEVICE Matrix4() {};
    AMREX_GPU_HOST_DEVICE Matrix4(Set::Scalar a_lambda, Set::Scalar a_mu) : lambda(a_lambda), mu(a_mu) {};
    AMREX_FORCE_INLINE
    Scalar operator () (const int i, const int j, const int k, const int l)
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
    void Print (std::ostream& /* os */)
    {
        Util::Message(INFO,"lambda = ", lambda, " mu = ", mu);
    }
    static Matrix4<3,Sym::Isotropic> Zero()
    {
        Matrix4<3,Sym::Isotropic> zero;
        zero.lambda = 0.0;
        zero.mu = 0.0;
        return zero;
    }
    friend Set::Matrix operator * (Matrix4<3,Sym::Isotropic> a, Set::Matrix b);
    AMREX_GPU_HOST_DEVICE void operator =  (Matrix4<3,Sym::Isotropic> a) {lambda =  a.lambda; mu =  a.mu;}
    AMREX_GPU_HOST_DEVICE void operator += (Matrix4<3,Sym::Isotropic> a) {lambda += a.lambda; mu += a.mu;}
    AMREX_GPU_HOST_DEVICE void operator -= (Matrix4<3,Sym::Isotropic> a) {lambda -= a.lambda; mu -= a.mu;}
    AMREX_GPU_HOST_DEVICE void operator *= (Matrix4<3,Sym::Isotropic> a) {lambda *= a.lambda; mu *= a.mu;}
    AMREX_GPU_HOST_DEVICE void operator /= (Matrix4<3,Sym::Isotropic> a) {lambda /= a.lambda; mu /= a.mu;}
    AMREX_GPU_HOST_DEVICE void operator *= (Set::Scalar alpha) {lambda *= alpha; mu *= alpha;}
    AMREX_GPU_HOST_DEVICE void operator /= (Set::Scalar alpha) {lambda /= alpha; mu /= alpha;}
};
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE Set::Matrix operator * (Matrix4<3,Sym::Isotropic> a, Set::Matrix b)
{
	Set::Matrix ret;
		
	ret(0,0) = (a.lambda + 2.*a.mu) * b(0,0) +       a.lambda      *b(1,1) +       a.lambda      *b(2,2);
	ret(1,1) =        a.lambda      * b(0,0) + (a.lambda + 2.*a.mu)*b(1,1) +       a.lambda      *b(2,2);
	ret(2,2) =        a.lambda      * b(0,0) +       a.lambda      *b(1,1) + (a.lambda + 2.*a.mu)*b(2,2);

	ret(1,2) = a.mu*(b(1,2) + b(2,1)); ret(2,1) = ret(1,2);
	ret(2,0) = a.mu*(b(2,0) + b(0,2)); ret(0,2) = ret(2,0);
	ret(0,1) = a.mu*(b(0,1) + b(1,0)); ret(1,0) = ret(0,1);
		
	return ret;
}
}
#endif