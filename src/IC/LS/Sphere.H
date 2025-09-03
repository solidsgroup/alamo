#ifndef IC_LS_SPHERE_H
#define IC_LS_SPHERE_H

#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

/// Initialize a spherical inclusion
namespace IC 
{
namespace LS 
{
class Sphere : public IC<Set::Scalar> {
public:
    static constexpr const char* name = "sphere";

    enum Type 
    { 
        XYZ, 
        YZ, 
        ZX, 
        XY 
    };

    // Constructors
    Sphere(amrex::Vector<amrex::Geometry>& _geom): IC(_geom) {};
    Sphere(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp): IC(_geom)
    {
        pp_queryclass(*this);
    };
    Sphere(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name): IC(_geom)
    {
        pp_queryclass(name, *this);
    };
    /// Constructor defining radius, and center coords
    Sphere(amrex::Vector<amrex::Geometry>& _geom, Set::Scalar _radius, Set::Vector _center, Type _type = Type::XYZ)
        : IC(_geom)
    {
        Define(_radius, _center, _type);
    }

    void Define(Set::Scalar a_radius,
        Set::Vector a_center,
        Type a_type)
    {
        radius = a_radius;
        center = a_center;
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

                Set::Scalar distance = std::sqrt(rsq) - radius;
                field(i, j, k, 0) = std::clamp(distance, InnerTube, OuterTube);
            });
        }
    };
    // This is a somewhat antiquated IC that will eventually be replaced
    // with the Expression IC.
    static void Parse(Sphere& value, IO::ParmParse& pp)
    {
        pp_query_default("radius", value.radius, 1.0); // Radius of the sphere
        pp_queryarr("center", value.center); // Vector location of the sphere center
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
    Type type = Type::XYZ;
    Set::Scalar TUBEWIDTH = 6.0;
};

} // namespace LS
} // namespace IC

#endif // IC_SPHERELS_H
