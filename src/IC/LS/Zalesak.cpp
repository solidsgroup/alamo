#include "IC/LS/Zalesak.H"

#include <cmath>
#include "IC/IC.H"
#include "Util/Util.H"

namespace IC {
namespace LS {

Zalesak::Zalesak(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

Zalesak::Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp) : IC(_geom) {
    pp_queryclass(*this);
}

Zalesak::Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, const std::string& name) : IC(_geom) {
    pp_queryclass(name, *this);
}

void Zalesak::Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar) {
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
            Set::Scalar slot_ls = computeSlotLevelSet(coords);
            Set::Scalar combined_distance = std::max(distance, -slot_ls); // Ensure subtraction
            field(i, j, k, 0) = computeLevelSetValue(combined_distance, Narrow_Band_Width, InnerTube, OuterTube);
        });
    }
}

void Zalesak::Parse(Zalesak& value, IO::ParmParse& pp) {
    pp_query("radius", value.radius);
    pp_queryarr("center", value.center);
    pp_queryarr("slot_center", value.slot_center);
    pp_queryarr("slot_size", value.slot_size);
    pp_query("fillet_radius", value.fillet_radius);
    std::string type_str;
    pp_query("type", type_str);

    if (type_str == "yz" || type_str == "zy") value.type = Type::YZ;
    else if (type_str == "zx" || type_str == "xz") value.type = Type::ZX;
    else if (type_str == "xy" || type_str == "yx") value.type = Type::XY;
    else value.type = Type::XYZ;
}

// Compute physical coordinates of the grid cell
AMREX_GPU_HOST_DEVICE Set::Vector Zalesak::computeCoordinates(int i, int j, int k, int lev, bool cellcentered) const {
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

Set::Scalar Zalesak::computeRSquared(const Set::Vector& coords, Type type, const Set::Vector& center) const {
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

Set::Scalar Zalesak::computeSlotLevelSet(const Set::Vector& coords) const {
    Set::Vector d = (coords - slot_center).cwiseAbs() - slot_size / 2.0;
    Set::Vector d_clamped = d.cwiseMax(Set::Vector::Zero());
    Set::Scalar outside_distance = d_clamped.norm();
    Set::Scalar inside_distance = std::min(std::max(d(0), std::max(d(1), d(2))), 0.0);
    return outside_distance + inside_distance;
}

/*Set::Scalar Zalesak::computeSlotLevelSet(const Set::Vector& coords) const {
    // Compute distance from rectangle edges
    Set::Vector d = (coords - slot_center).cwiseAbs() - slot_size / 2.0;
    
    // If inside the rectangle, return negative distance
    if (d(0) <= 0 && d(1) <= 0) {
        return -std::sqrt(d.squaredNorm());
    }

    // Handle fillet at the corners
    Set::Scalar corner_dist = std::sqrt(std::max(d(0), 0.0) * std::max(d(0), 0.0) +
                                        std::max(d(1), 0.0) * std::max(d(1), 0.0));

    // Apply fillet smoothing
    //std::cout << "fillet_radius:" << " " << fillet_radius << std::endl;
    return (corner_dist - fillet_radius);
}*/

Set::Scalar Zalesak::computeLevelSetValue(Set::Scalar distance, Set::Scalar Narrow_Band_Width, Set::Scalar InnerTube, Set::Scalar OuterTube) const {
    if (std::abs(distance) <= Narrow_Band_Width) return distance;
    if (distance < -Narrow_Band_Width) return InnerTube;
    if (distance > Narrow_Band_Width) return OuterTube;
    return NAN;
}

} // namespace LS
} // namespace IC
