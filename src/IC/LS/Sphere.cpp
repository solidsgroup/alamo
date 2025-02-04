#include "IC/LS/Sphere.H"

#include <cmath>
#include "IC/IC.H"
#include "Util/Util.H"

namespace IC {
namespace LS {

// Constructors
Sphere::Sphere(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

Sphere::Sphere(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp) : IC(_geom) {
    pp_queryclass(*this);
}

Sphere::Sphere(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, const std::string& name) : IC(_geom) {
    pp_queryclass(name, *this);
}

Sphere::Sphere(amrex::Vector<amrex::Geometry>& _geom, Set::Scalar _radius, Set::Vector _center, Type _type)
    : IC(_geom), radius(_radius), center(_center), type(_type) {}

// Add the level set field
void Sphere::Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar) {
    bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));
    const Set::Scalar* DX = geom[lev].CellSize();
    const Set::Scalar Narrow_Band_Width = 6.0 * DX[0];
    const Set::Scalar InnerTube = -Narrow_Band_Width;
    const Set::Scalar OuterTube = Narrow_Band_Width;

    for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox();
        auto field = a_field[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto coords = computeCoordinates(i, j, k, lev, cellcentered);
            Set::Scalar rsq = computeRSquared(coords, type, center);
            Set::Scalar distance = std::sqrt(rsq) - radius;
            field(i, j, k, 0) = computeLevelSetValue(distance, Narrow_Band_Width, InnerTube, OuterTube);
        });
    }
}

// Parse function to set parameters from input
void Sphere::Parse(Sphere& value, IO::ParmParse& pp) {
    pp_query("radius", value.radius);
    pp_queryarr("center", value.center);
    std::string type_str;
    pp_query("type", type_str);

    if (type_str == "yz" || type_str == "zy") value.type = Type::YZ;
    else if (type_str == "zx" || type_str == "xz") value.type = Type::ZX;
    else if (type_str == "xy" || type_str == "yx") value.type = Type::XY;
    else value.type = Type::XYZ;
}

// Compute physical coordinates of the grid cell
AMREX_GPU_HOST_DEVICE inline Set::Vector Sphere::computeCoordinates(int i, int j, int k, int lev, bool cellcentered) const {
    Set::Vector coords;
    const auto& DX = geom[lev].CellSize();
    const auto& prob_lo = geom[lev].ProbLo();

    if (cellcentered) {
        AMREX_D_TERM(coords(0) = prob_lo[0] + (i + 0.5) * DX[0];,
                     coords(1) = prob_lo[1] + (j + 0.5) * DX[1];,
                     coords(2) = prob_lo[2] + (k + 0.5) * DX[2];);
    } else {
        AMREX_D_TERM(coords(0) = prob_lo[0] + i * DX[0];,
                     coords(1) = prob_lo[1] + j * DX[1];,
                     coords(2) = prob_lo[2] + k * DX[2];);
    }
    return coords;
}

// Compute rsq based on type
Set::Scalar Sphere::computeRSquared(const Set::Vector& coords, Type type, const Set::Vector& center) const {
    switch (type) {
    case Type::XYZ:
        return AMREX_D_TERM((coords(0) - center(0)) * (coords(0) - center(0)),
                            + (coords(1) - center(1)) * (coords(1) - center(1)),
                            + (coords(2) - center(2)) * (coords(2) - center(2)));

    case Type::XY:
        return (coords(0) - center(0)) * (coords(0) - center(0)) +
               (coords(1) - center(1)) * (coords(1) - center(1));

    case Type::YZ:
        return (coords(1) - center(1)) * (coords(1) - center(1)) +
               (coords(2) - center(2)) * (coords(2) - center(2));

    case Type::ZX:
        return (coords(0) - center(0)) * (coords(0) - center(0)) +
               (coords(2) - center(2)) * (coords(2) - center(2));

    default:
        return NAN;
    }
}

// Compute level set value
AMREX_GPU_HOST_DEVICE inline Set::Scalar Sphere::computeLevelSetValue(Set::Scalar distance, Set::Scalar Narrow_Band_Width, Set::Scalar InnerTube, Set::Scalar OuterTube) const {
    if (std::abs(distance) <= Narrow_Band_Width) return distance;
    if (distance < -Narrow_Band_Width) return InnerTube;
    if (distance > Narrow_Band_Width) return OuterTube;
    return NAN;
}

} // namespace LS
} // namespace IC
