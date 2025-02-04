#include "IC/LS/Zalesak.H"

#include <cmath>
#include "IC/IC.H"
#include "Util/Util.H"

namespace IC {
namespace LS {

void Zalesak::Define(Set::Scalar a_radius, Set::Vector a_center, Set::Scalar a_slot_width, Set::Scalar a_slot_length, Type a_type) {
    radius = a_radius;
    center = a_center;
    type = a_type;
    slot_width = (a_slot_width > 0) ? a_slot_width : 0.1 * radius;
    slot_length = (a_slot_length > 0) ? a_slot_length : 0.6 * radius;
}

void Zalesak::Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar) {
    bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));
    const Set::Scalar* DX = geom[lev].CellSize();
    
    for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.tilebox();
        amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
        
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar x, y, z;
            if (cellcentered) {
                x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * DX[0];
                y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * DX[1];
                z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * DX[2];
            } else {
                x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * DX[0];
                y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * DX[1];
                z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * DX[2];
            }

            Set::Scalar rsq = (type == Type::XYZ) ? (x - center(0)) * (x - center(0)) + (y - center(1)) * (y - center(1)) + (z - center(2)) * (z - center(2))
                         : (type == Type::XY)  ? (x - center(0)) * (x - center(0)) + (y - center(1)) * (y - center(1))
                         : (type == Type::YZ)  ? (y - center(1)) * (y - center(1)) + (z - center(2)) * (z - center(2))
                         : (x - center(0)) * (x - center(0)) + (z - center(2)) * (z - center(2));
            
            Set::Scalar distance = std::sqrt(rsq) - radius;
            if (x > center(0) - 0.5 * slot_width && x < center(0) + 0.5 * slot_width &&
                y < center(1) + radius && y > center(1) + radius - slot_length) {
                field(i, j, k, 0) = 1.0;
            } else {
                field(i, j, k, 0) = distance;
            }
        });
    }
}

void Zalesak::Parse(Zalesak& value, IO::ParmParse& pp) {
    pp_query("radius", value.radius);
    pp_queryarr("center", value.center);
    pp_query("slot_width", value.slot_width);
    pp_query("slot_length", value.slot_length);
    std::string type;
    pp_query("type", type);
    if (type == "yz" || type == "zy")      value.type = Type::YZ;
    else if (type == "zx" || type == "xz") value.type = Type::ZX;
    else if (type == "xy" || type == "yx") value.type = Type::XY;
    else                                   value.type = Type::XYZ;
}

} // namespace LS
} // namespace IC
