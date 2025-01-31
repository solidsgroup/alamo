#ifndef IC_LS_ZALESAKLS_H_
#define IC_LS_ZALESAKLS_H_

#include <cmath>
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace IC {
namespace LS {

class Zalesak : public IC {
public:
    enum Type {
        XYZ,
        YZ,
        ZX,
        XY
    };

    Zalesak(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {} ;
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp) : IC(_geom)
    {pp_queryclass(*this);}
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {pp_queryclass(name,*this);}
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, Set::Scalar _radius, Set::Vector _center,
              Set::Scalar _slot_width = -1, Set::Scalar _slot_length = -1, Type _type = Type::XYZ);

    void Define(Set::Scalar a_radius, Set::Vector a_center, Set::Scalar a_slot_width = -1, Set::Scalar a_slot_length = -1, Type a_type = Type::XYZ);
    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar);
    static void Parse(Zalesak& value, IO::ParmParse& pp);

private:
    Set::Vector center = Set::Vector::Zero();
    Set::Scalar radius = 1.0;
    Set::Scalar slot_width = 0.1;
    Set::Scalar slot_length = 0.6;
    Type type = Type::XYZ;
};

} // namespace LS
} // namespace IC

#endif // IC_LS_ZALESAKLS_H_
