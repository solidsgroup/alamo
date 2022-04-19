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
    void Copy(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) {}
    int NComp() {return 1;}
};

template<>
inline void Field<Set::Vector>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost)
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
inline int Field<Set::Vector>::NComp() {return AMREX_SPACEDIM;} 

template<>
inline void Field<Set::Matrix>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost)
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
inline int Field<Set::Matrix>::NComp() {return AMREX_SPACEDIM*AMREX_SPACEDIM;} 



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

    void Copy(int /*a_lev*/, amrex::MultiFab &/*a_dst*/, int /*a_dstcomp*/, int /*a_nghost*/) 
    {Util::Abort(INFO,"This should never get called");}
    int NComp() 
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
