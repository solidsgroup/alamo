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
{};

template <>
class Field<Set::Scalar> : public amrex::Vector<std::unique_ptr<amrex::MultiFab>>
{
public:
    void Define(int a_levs, amrex::Vector<amrex::BoxArray> a_grids, amrex::Vector<amrex::DistributionMapping> a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_levs == a_grids.size()));
        Util::Assert(INFO,TEST(a_levs == a_dmap.size()));
        this->resize(a_levs);
        for (int lev = 0; lev < a_levs; lev++)
        {
            (*this)[lev].reset(new amrex::MultiFab(a_grids[lev],a_dmap[lev],a_ncomp,a_nghost));
        }
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
