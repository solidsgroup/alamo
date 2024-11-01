#ifndef IC_SPHERELS_H_
#define IC_SPHERELS_H_

#include <cmath>

#include "IC/IC.H"
#include "Util/Util.H"

namespace IC
{
class SphereLS: public IC
{
public:
    enum Type
    {
        XYZ,
        YZ,
        ZX,
        XY
    };

    SphereLS(amrex::Vector<amrex::Geometry>& _geom): IC(_geom) {}
    SphereLS(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp): IC(_geom)
    {
        pp_queryclass(*this);
    }
    SphereLS(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name): IC(_geom)
    {
        pp_queryclass(name, *this);
    }
    SphereLS(amrex::Vector<amrex::Geometry>& _geom, Set::Scalar _radius, Set::Vector _center, Type _type = Type::XYZ)
        : IC(_geom)
    {
        Define(_radius, _center, _type);
    }

    void Define(Set::Scalar a_radius, Set::Vector a_center, Type a_type)
    {
        radius = a_radius;
        center = a_center;
        type = a_type;
    }

    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar)
    {
        bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));
        int ncomp = a_field[0]->nComp();

        // Get the cell size corresponding to this level
        const Set::Scalar* DX = geom[lev].CellSize();
        /*amrex::Print() << "DX is" << DX[0] << "\n";*/

        // Narrow band width for level set
        const Set::Scalar Narrow_Band_Width = 6.0 * DX[0]; // Assuming DX[0] for simplicity in defining narrow_band_width
        const Set::Scalar InnerTube = -Narrow_Band_Width;
        const Set::Scalar OuterTube = Narrow_Band_Width;

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar AMREX_D_DECL(x, y, z);
                if (cellcentered)
                {
                    AMREX_D_TERM(
                        x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                        y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                        z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];
                    );
                }
                else
                {
                    AMREX_D_TERM(
                        x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * geom[lev].CellSize()[0];,
                        y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * geom[lev].CellSize()[1];,
                        z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * geom[lev].CellSize()[2];
                    );
                }

                Set::Scalar rsq = NAN;

                if (type == Type::XYZ)
                {
                    // 3D Sphere
                    rsq = AMREX_D_TERM(
                        (x - center(0)) * (x - center(0)),
                        +(y - center(1)) * (y - center(1)),
                        +(z - center(2)) * (z - center(2))
                    );
                }
                else if (type == Type::XY)
                {
                    // Cylinder along Z axis
                    rsq = AMREX_D_TERM(
                        (x - center(0)) * (x - center(0)),
                        +(y - center(1)) * (y - center(1)),
                        
                    );
                }
                else if (type == Type::YZ)
                {
                    // Cylinder along X axis
                    rsq = AMREX_D_TERM(
                        ,
                        +(y - center(1)) * (y - center(1)),
                        +(z - center(2)) * (z - center(2))
                    );
                }
                else if (type == Type::ZX)
                {
                    // Cylinder along Y axis
                    rsq = AMREX_D_TERM(
                        (x - center(0)) * (x - center(0)),
                        ,
                        +(z - center(2)) * (z - center(2))
                    );
                }

                Set::Scalar distance = std::sqrt(rsq) - radius;

                if (std::abs(distance) <= Narrow_Band_Width)
                {
                    // Initialize narrow band level set
                    field(i, j, k, 0) = distance;
                    /*amrex::Print() << "field(" << i << ", " << j << ", " << k << ") = " << field(i, j, k) << "\n";*/
                }
                else
                {
                    // Inside the interface and beyond the narrow band (INNERTUBE)
                    if (distance < -Narrow_Band_Width)
                    {
                        field(i, j, k, 0) = InnerTube; // Set to -6 * DX for INNERTUBE
                    }
                    // Outside the interface and beyond the narrow band (OUTERTUBE)
                    else if (distance > Narrow_Band_Width)
                    {
                        field(i, j, k, 0) = OuterTube; // Set to +6 * DX for OUTERTUBE
                    }
                }
            });
        }
    }

    static void Parse(SphereLS& value, IO::ParmParse& pp)
    {
        pp_query("radius", value.radius);
        pp_queryarr("center", value.center);
        std::string type;
        pp_query("type", type);
        if (type == "yz")
            value.type = Type::YZ;
        if (type == "zx")
            value.type = Type::ZX;
        if (type == "xy")
            value.type = Type::XY;
        if (type == "xyz")
            value.type = Type::XYZ;
    }

private:
    Set::Vector center = Set::Vector::Zero();
    Set::Scalar radius = 1.0;
    Type type = Type::XYZ;
};
}

#endif // IC_SPHERELS_H_

