#ifndef SET_MATRIX3_H
#define SET_MATRIX3_H

#include "Set/Base.H"

namespace Set
{
class Matrix3
{
public:
    Set::Matrix data[AMREX_SPACEDIM];
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Scalar & operator () (const int i, const int j, const int k)
    {
        return data[i](j,k);
    }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    const Scalar & operator () (const int i, const int j, const int k) const
    {
        return data[i](j,k);
    }
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Matrix & operator [] (const int i)
    {
        return data[i];
    }

    //AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void operator  = (Matrix3 a) {for (int i = 0; i < AMREX_SPACEDIM; i++) data[i]  = a.data[i];}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void operator += (Matrix3 a) {for (int i = 0; i < AMREX_SPACEDIM; i++) data[i] += a.data[i];}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void operator -= (Matrix3 a) {for (int i = 0; i < AMREX_SPACEDIM; i++) data[i] -= a.data[i];}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void operator *= (Set::Scalar alpha) {for (int i = 0; i < AMREX_SPACEDIM; i++) data[i] *= alpha;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void operator /= (Set::Scalar alpha) {for (int i = 0; i < AMREX_SPACEDIM; i++) data[i] /= alpha;}

    static Matrix3 Zero()
    {
        Matrix3 ret;
        for (int i = 0; i < AMREX_SPACEDIM; i++) ret.data[i] = Set::Matrix::Zero();
        return ret;
    }
    static Matrix3 Random()
    {
        Matrix3 ret;
        for (int i = 0; i < AMREX_SPACEDIM; i++) ret.data[i] = Set::Matrix::Random();
        return ret;
    }

};
}
#endif
