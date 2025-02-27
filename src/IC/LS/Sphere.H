#ifndef IC_LS_SPHERE_H
#define IC_LS_SPHERE_H

#include <cmath>
#include "IC/IC.H"
#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace IC {
namespace LS {

class Sphere : public IC<Set::Scalar> {
public:
    static constexpr const char* name = "sphere";

    enum Type { XYZ, YZ, ZX, XY };

    // Constructors
    explicit Sphere(amrex::Vector<amrex::Geometry>& _geom);
    Sphere(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp);
    Sphere(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, const std::string& name);
    Sphere(amrex::Vector<amrex::Geometry>& _geom, Set::Scalar _radius, Set::Vector _center, Type _type = Type::XYZ);

    // Add the level set field
    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar);

    // Parse function to set parameters from input
    static void Parse(Sphere& value, IO::ParmParse& pp);

private:
    Set::Scalar radius = 1.0;
    Set::Vector center = Set::Vector::Zero();
    Type type = Type::XYZ;

    // Compute physical coordinates of the grid cell
    AMREX_GPU_HOST_DEVICE inline Set::Vector computeCoordinates(int i, int j, int k, int lev, bool cellcentered) const;

    // Compute rsq based on type
    AMREX_GPU_HOST_DEVICE inline Set::Scalar computeRSquared(const Set::Vector& coords, Type type, const Set::Vector& center) const;

    // Compute level set value
    AMREX_GPU_HOST_DEVICE inline Set::Scalar computeLevelSetValue(Set::Scalar distance, Set::Scalar Narrow_Band_Width, Set::Scalar InnerTube, Set::Scalar OuterTube) const;
};

} // namespace LS
} // namespace IC

#endif // IC_SPHERELS_H
