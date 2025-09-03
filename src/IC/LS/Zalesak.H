#ifndef IC_LS_ZALESAK_H
#define IC_LS_ZALESAK_H

#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

/// Initialize a zalesak disk
namespace IC 
{
namespace LS 
{
class Zalesak : public IC<Set::Scalar> {
public:
    static constexpr const char* name = "zalesak";

    enum Type 
    { 
        XYZ, 
        YZ, 
        ZX, 
        XY 
    };

    //  Constructors
    Zalesak(amrex::Vector<amrex::Geometry>& _geom): IC(_geom) {};
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp): IC(_geom)
    {
        pp_queryclass(*this);
    }
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, const std::string& name): IC(_geom)
    {
        pp_queryclass(name, *this);
    }
    /// Constructor defining radius, and center coords
    Zalesak(amrex::Vector<amrex::Geometry>& _geom, Set::Scalar _radius, Set::Vector _center, Set::Vector _slot_center, Set::Vector _slot_size,  Set::Scalar _fillet_radius, Type _type = Type::XYZ)
        : IC(_geom)
    {
        Define(_radius, _center, _slot_center, _slot_size, _fillet_radius, _type);
    }

    void Define(Set::Scalar a_radius,
        Set::Vector a_center,
        Set::Vector a_slot_center,
        Set::Vector a_slot_size,
        Set::Scalar a_fillet_radius,
        Type a_type)
    {
        radius = a_radius;
        center = a_center;
        slot_center = a_slot_center,
        slot_size = a_slot_size,
        fillet_radius = a_fillet_radius,
        type = a_type;
    };

    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar)
    {
        bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));
        int nghost = a_field[0]->nGrow();

        const Set::Scalar* DX = geom[lev].CellSize();
        const Set::Scalar min_DX = *std::min_element(DX, DX + AMREX_SPACEDIM);
        const Set::Scalar Narrow_Band_Width = TUBEWIDTH * min_DX;
        const Set::Scalar InnerTube = -Narrow_Band_Width;
        const Set::Scalar OuterTube = Narrow_Band_Width;

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::IndexType index_type = a_field[lev]->ixType();
            if (index_type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox(nghost);
            else if (index_type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox(nghost);
            else Util::Abort(INFO, "Unknown index type");

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {

                Set::Scalar AMREX_D_DECL(x, y, z);
                if (cellcentered)
                {
                    AMREX_D_TERM(x = geom[lev].ProbLo()[0] + ((amrex::Real)(i)+0.5) * geom[lev].CellSize()[0];,
                        y = geom[lev].ProbLo()[1] + ((amrex::Real)(j)+0.5) * geom[lev].CellSize()[1];,
                        z = geom[lev].ProbLo()[2] + ((amrex::Real)(k)+0.5) * geom[lev].CellSize()[2];);
                }
                else
                {
                    AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i)*geom[lev].CellSize()[0];,
                        y = geom[lev].ProbLo()[1] + (amrex::Real)(j)*geom[lev].CellSize()[1];,
                        z = geom[lev].ProbLo()[2] + (amrex::Real)(k)*geom[lev].CellSize()[2];);
                }

                // Compute r^2 term for circular distance
                Set::Scalar rsq = NAN;

                if (type == Type::XYZ)
                {
                    // 3D Sphere
                    rsq =
                        AMREX_D_TERM((x - center(0)) * (x - center(0)),
                            +(y - center(1)) * (y - center(1)),
                            +(z - center(2)) * (z - center(2)));
                }
                else if (type == Type::XY)
                {
                    // Cylinder along Z axis
                    rsq =
                        AMREX_D_TERM((x - center(0)) * (x - center(0)),
                            +(y - center(1)) * (y - center(1)),
                            );
                }
                else if (type == Type::YZ)
                {
                    // Cylinder along X axis
                    rsq =
                        AMREX_D_TERM(,
                            +(y - center(1)) * (y - center(1)),
                            +(z - center(2)) * (z - center(2)));
                }
                else if (type == Type::ZX)
                {
                    // Cylinder along Y axis
                    rsq =
                        AMREX_D_TERM((x - center(0)) * (x - center(0)),
                            ,
                            +(z - center(2)) * (z - center(2)));
                }

                Set::Scalar circular_distance = std::sqrt(rsq) - radius;

                // Compute slot signed distance (rounded box SDF)
                Set::Vector d = Set::Vector::Zero();
                Set::Vector d_clamped = Set::Vector::Zero();

                d[0] = std::abs(x - slot_center(0)) - slot_size(0) / 2.0;
                #if (AMREX_SPACEDIM >= 2)
                d[1] = std::abs(y - slot_center(1)) - slot_size(1) / 2.0;
                #endif
                #if (AMREX_SPACEDIM == 3)
                d[2] = std::abs(z - slot_center(2)) - slot_size(2) / 2.0;
                #endif

                d_clamped[0] = std::max(d[0], 0.0);
                #if (AMREX_SPACEDIM >= 2)
                d_clamped[1] = std::max(d[1], 0.0);
                #endif
                #if (AMREX_SPACEDIM == 3)
                d_clamped[2] = std::max(d[2], 0.0);
                #endif

                Set::Scalar outside_dist = 0.0;
                outside_dist += d_clamped[0] * d_clamped[0];
                #if (AMREX_SPACEDIM >= 2)
                outside_dist += d_clamped[1] * d_clamped[1];
                #endif
                #if (AMREX_SPACEDIM == 3)
                outside_dist += d_clamped[2] * d_clamped[2];
                #endif

                outside_dist = std::sqrt(outside_dist);

                Set::Scalar inside_offset = d[0];
                #if (AMREX_SPACEDIM >= 2)
                inside_offset = std::max(inside_offset, d[1]);
                #endif
                #if (AMREX_SPACEDIM == 3)
                inside_offset = std::max(inside_offset, d[2]);
                #endif
                inside_offset = std::min(inside_offset, 0.0);

                Set::Scalar slot_distance = outside_dist + inside_offset - fillet_radius;

                Set::Scalar levelset = std::max(circular_distance, -slot_distance);
                field(i, j, k, 0) = std::clamp(levelset, InnerTube, OuterTube);
            });
        }
    };

    static void Parse(Zalesak& value, IO::ParmParse& pp)
    {
        pp_query_default("radius", value.radius, 1.0); // Radius of the sphere
        pp_queryarr("center", value.center); // Vector location of the sphere center
        pp_queryarr("slot_center", value.slot_center); // the center coordinate of the slot rectangle that is cut out
        pp_queryarr("slot_size", value.slot_size); // lengths in x, y, and z planes of rectangle that is cut out
        pp_query("fillet_radius", value.fillet_radius); // the radius of the fillet used to smooth slot
        std::string type;
        // Type - can be cylinder oriented along the x, y, z directions or full sphere.
        pp_query_validate("type", type, {"xyz","yz","zx","xy"}); 
        if (type == "yz")  value.type = Type::YZ;
        if (type == "zx")  value.type = Type::ZX;
        if (type == "xy")  value.type = Type::XY;
        if (type == "xyz") value.type = Type::XYZ;
    }

private:
    Set::Scalar radius = NAN;
    Set::Vector center = Set::Vector::Zero();
    Set::Vector slot_center = Set::Vector::Zero();
    Set::Vector slot_size = Set::Vector::Zero();
    Set::Scalar fillet_radius = NAN;
    Type type = Type::XYZ;
    Set::Scalar TUBEWIDTH = 6.0;
};

} // namespace LS
} // namespace IC

#endif // IC_LS_ZALESAK_H
