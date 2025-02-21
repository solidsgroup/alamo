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
/// Pure abstract IC object from which all other IC objects inherit.
template<class T = Set::Scalar>
class IC
{
public:
    IC (amrex::Vector<amrex::Geometry> &_geom)
        : geom(_geom) {} ;
    virtual ~IC() {}

    virtual void Add(const int &lev, Set::Field<T> &field, Set::Scalar time) = 0;
    void Add(const int &lev, Set::Field<T> &field)
    {
        Add(lev,field,0.0);
    }


    // SFINAE solution: this is how we initialize if the type is arithmetic and has a  
    // "0" value
    template <typename U = T, typename std::enable_if_t<std::is_arithmetic_v<U>, int> = 0>
    void Initialize(const int &a_lev,
                    Set::Field<T> &a_field,
                    Set::Scalar a_time = 0.0)
    {
        Util::Assert(INFO,TEST(a_lev < a_field.size())," a_lev=",a_lev," size=",a_field.size());
        a_field[a_lev]->setVal(0.0);
        Add(a_lev,a_field,a_time);
    };

    // SFINAE solution: we initialize using "Zero" if it is a class type.
    template <typename U = T, typename std::enable_if_t<!std::is_arithmetic_v<U>, int> = 0>
    void Initialize(const int &a_lev,
                    Set::Field<T> &a_field,
                    Set::Scalar a_time = 0.0)
    {
        Util::Assert(INFO,TEST(a_lev < a_field.size())," a_lev=",a_lev," size=",a_field.size());
        a_field[a_lev]->setVal(T::Zero());
        Add(a_lev,a_field,a_time);
    };

    virtual void SetComp(int a_comp) final {comp = a_comp;}
protected:
    amrex::Vector<amrex::Geometry> &geom;
    int comp = 0;
};
}
#endif
