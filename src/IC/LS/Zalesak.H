#ifndef IC_LS_ZALESAK_H
#define IC_LS_ZALESAK_H

#include <cmath>
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace IC {
namespace LS {

class Zalesak : public IC {
public:
    static constexpr const char* name = "zalesak";

    enum Type { XYZ, YZ, ZX, XY };

    explicit Zalesak(amrex::Vector<amrex::Geometry>& _geom);
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp);
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, const std::string& name);

    static void Parse(Zalesak& value, IO::ParmParse& pp);
    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar);
    
private:
    Set::Scalar radius = 1.0;
    Set::Vector center = Set::Vector::Zero();
    Set::Vector slot_center = Set::Vector::Zero();
    Set::Vector slot_size = Set::Vector::Zero();
    Set::Scalar fillet_radius = 0.0;
    Type type = Type::XYZ;

    AMREX_GPU_HOST_DEVICE inline Set::Vector computeCoordinates(int i, int j, int k, int lev, bool cellcentered) const;
    AMREX_GPU_HOST_DEVICE inline Set::Scalar computeCircleLevelSet(const Set::Vector& coords) const;
    AMREX_GPU_HOST_DEVICE inline Set::Scalar computeRSquared(const Set::Vector& coords, Type type, const Set::Vector& center) const;
    AMREX_GPU_HOST_DEVICE inline Set::Scalar computeSlotLevelSet(const Set::Vector& coords) const;
    AMREX_GPU_HOST_DEVICE inline Set::Scalar computeLevelSetValue(Set::Scalar distance, Set::Scalar Narrow_Band_Width, Set::Scalar InnerTube, Set::Scalar OuterTube) const;
};

} // namespace LS
} // namespace IC

#endif // IC_LS_ZALESAK_H

