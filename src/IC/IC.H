//
// Initial condition (IC) objects are used to set fields to specified values,
// through mathematical expressions, images, etc.
// They are often used as initial conditions, but may be used generally to
// incorporate static or time-dependent values.
// All IC methods inherit from the base IC class, and are implemented by
// overriding the "Add" method.
// (Initialize uses Add after setting the field to zero.)
//
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

    virtual void Add(const int &lev, Set::Field<Set::Scalar> &field, Set::Scalar time) = 0;
    virtual void Add(const int &lev, Set::Field<Set::Scalar> &field)
    {
        Add(lev,field,0.0);
    }
    virtual void Add(const int &, Set::Field<Set::Vector> &, Set::Scalar) 
    {Util::Abort(INFO,"Not yet implemented");};
    void Initialize(const int &a_lev,
                    Set::Field<Set::Scalar> &a_field,
                    Set::Scalar a_time = 0.0)
    {
        Util::Assert(INFO,TEST(a_lev < a_field.size())," a_lev=",a_lev," size=",a_field.size());
        a_field[a_lev]->setVal(0.0);
        Add(a_lev,a_field,a_time);
    };
    void Initialize(const int &a_lev,
                    Set::Field<Set::Vector> &a_field,
                    Set::Scalar a_time = 0.0)
    {
        Util::Assert(INFO,TEST(a_lev < a_field.size())," a_lev=",a_lev," size=",a_field.size());
        a_field[a_lev]->setVal(Set::Vector::Zero());
        Add(a_lev,a_field,a_time);
    };

    virtual void SetComp(int a_comp) final {comp = a_comp;}
protected:
    amrex::Vector<amrex::Geometry> &geom;
    int comp = 0;
};
}
#endif
