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
};

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
};
}

namespace Util
{
Set::Scalar Random();
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
