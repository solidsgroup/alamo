#ifndef IC_IC_H_
#define IC_IC_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Util/Util.H"
#include "Set/Set.H"

namespace IC
{
/// \brief Pure abstract IC object from which all other IC objects inherit.
class IC
{
public:
    IC (amrex::Vector<amrex::Geometry> &_geom)
        : geom(_geom) {} ;
    virtual ~IC() {}

    virtual void Add(const int &lev, Set::Field<Set::Scalar> &field) = 0;
    void Initialize(const int &a_lev,
                    Set::Field<Set::Scalar> &a_field)
    {
        a_field[a_lev]->setVal(0.0);
        Add(a_lev,a_field);
    };

    virtual void SetComp(int a_comp) final {comp = a_comp;}
protected:
    amrex::Vector<amrex::Geometry> &geom;
    int comp = 0;
};
}
#endif
