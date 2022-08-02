#ifndef SET_SET_
#define SET_SET_

#include <iomanip>

#include "Util/Util.H"

#include "Set/Base.H"
/// \brief A collection of data types and symmetry-reduced data structures

namespace Set
{
//class Field : public amrex::MultiFab
//{

//};
template <class T>
class Field : public amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>>
{
public:
    Field() {} 
    Field(int a_levs, const amrex::Vector<amrex::BoxArray> & a_grids, const amrex::Vector<amrex::DistributionMapping> & a_dmap, int a_ncomp, int a_nghost) 
    {
        Define(a_levs,a_grids,a_dmap,a_ncomp,a_nghost);
    } 
    Field(int size) : amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>>(size) {}
    void Define(int a_levs, const amrex::Vector<amrex::BoxArray> & a_grids, const amrex::Vector<amrex::DistributionMapping> & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_levs == a_grids.size()));
        Util::Assert(INFO,TEST(a_levs == a_dmap.size()));
        this->resize(a_levs);
        for (int lev = 0; lev < a_levs; lev++)
        {
            (*this)[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(a_grids[lev],a_dmap[lev],a_ncomp,a_nghost));
        }
    }
    void Define(int a_lev, const amrex::BoxArray & a_grid, const amrex::DistributionMapping & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_lev < this->size()));
        (*this)[a_lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(a_grid,a_dmap,a_ncomp,a_nghost));
    }
    int finest_level = 0;
    void Copy(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) const {}
    void CopyFrom(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) const {}
    void Add(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) const {}
    void AddFrom(int /*a_lev*/, amrex::MultiFab &/*a_src*/, int /*a_srccomp*/, int /*a_nghost*/) const {}
    int NComp() const {return 0;}
    virtual std::string Name(int) const {return name;}
    std::string name;    
};

template<>
inline void Field<Set::Vector>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Set::Vector> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k)(n);
                });
            }
        }
    }    
}
template<>
inline void Field<Set::Vector>::CopyFrom(int a_lev, amrex::MultiFab &a_src, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<Set::Vector> const & dst = ((*this)[a_lev])->array(mfi);
            amrex::Array4<const Set::Scalar> const & src = a_src.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k)(n) = src(i,j,k,a_dstcomp + n);
                });
            }
        }
    }    
}
template<>
inline void Field<Set::Vector>::Add(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Set::Vector> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) += src(i,j,k)(n);
                });
            }
        }
    }    
}
template<>
inline void Field<Set::Vector>::AddFrom(int a_lev, amrex::MultiFab &a_src, int a_srccomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_src, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<Set::Vector> const & dst = ((*this)[a_lev])->array(mfi);
            amrex::Array4<const Set::Scalar> const & src = a_src.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k)(n) += src(i,j,k,a_srccomp + n);
                });
            }
        }
    }    
}

template<> inline int         Field<Set::Vector>::NComp() const {return AMREX_SPACEDIM;} 
template<> inline std::string Field<Set::Vector>::Name(int i) const
{
    #if AMREX_SPACEDIM>0
    if (i==0) return name + "_x";
    #endif
    #if AMREX_SPACEDIM>1
    else if (i==1) return name + "_y";
    #endif
    #if AMREX_SPACEDIM>2
    else if (i==2) return name + "_z";
    #endif
    else Util::Abort(INFO,"Invalid component");
    return "";
}

template<>
inline void Field<Set::Matrix>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Set::Matrix> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM*AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k)(n/AMREX_SPACEDIM,n%AMREX_SPACEDIM);
                });
            }
        }
    }    
}
template<>
inline int Field<Set::Matrix>::NComp() const {return AMREX_SPACEDIM*AMREX_SPACEDIM;} 
template<> inline std::string Field<Set::Matrix>::Name(int i) const
{
    #if AMREX_SPACEDIM==1
    if (i==0) return name + "_xx";
    #elif AMREX_SPACEDIM==2
    if (i==0)      return name + "_xx";
    else if (i==1) return name + "_xy";
    else if (i==2) return name + "_yx";
    else if (i==3) return name + "_yy";
    #elif AMREX_SPACEDIM==3
    if (i==0)      return name + "_xx";
    else if (i==1) return name + "_xy";
    else if (i==2) return name + "_xz";
    else if (i==3) return name + "_yx";
    else if (i==4) return name + "_yy";
    else if (i==5) return name + "_yz";
    else if (i==6) return name + "_zx";
    else if (i==7) return name + "_zy";
    else if (i==8) return name + "_zz";
    #endif
    else Util::Abort(INFO, "Invalid component");
    return "";
}


template <>
class Field<Set::Scalar> : public amrex::Vector<std::unique_ptr<amrex::MultiFab>>
{
public:
    Field() {} 
    Field(int size) : amrex::Vector<std::unique_ptr<amrex::MultiFab>>(size) {}
    void Define(int a_levs, const amrex::Vector<amrex::BoxArray> & a_grids, const amrex::Vector<amrex::DistributionMapping> & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_levs == a_grids.size()));
        Util::Assert(INFO,TEST(a_levs == a_dmap.size()));
        this->resize(a_levs);
        for (int lev = 0; lev < a_levs; lev++)
        {
            (*this)[lev].reset(new amrex::MultiFab(a_grids[lev],a_dmap[lev],a_ncomp,a_nghost));
        }
    }
    void Define(int a_lev, const amrex::BoxArray & a_grid, const amrex::DistributionMapping & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_lev < this->size()));
        (*this)[a_lev].reset(new amrex::MultiFab(a_grid,a_dmap,a_ncomp,a_nghost));
    }
    int finest_level = 0;

    void Copy(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) const
    {Util::Abort(INFO,"This should never get called");}
    int NComp() const
    {Util::Abort(INFO,"This should never be called");return -1;}
};

}

namespace Util
{
Set::Scalar Random();
Set::Scalar Gaussian(amrex::Real mean,amrex::Real std_deviation);
}

namespace Set
{
namespace Constant
{
static const Set::Scalar Pi = 3.14159265359;
}
}

#include "Set/Matrix3.H"
#include "Set/Matrix4.H"


#endif
