//
// Fill a domain (region where field=0) with packed spheres (regions where field=1).
// Sphere locations and radii are determined from an xyzr file.
//

#ifndef IC_PSREAD_H_
#define IC_PSREAD_H_

#include "Set/Set.H"
#include "IC/IC.H"
#include "Unit/Unit.H"
#include "AMReX_GpuAsyncArray.H"
#include "AMReX_ParallelDescriptor.H"
#include "Util/Util.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <mpi.h>
#include "IO/ParmParse.H"

namespace IC
{
class PSRead : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "psread";

    enum Type
    {
        Partition,
        Values
    };

    PSRead(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PSRead(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp.queryclass(name, *this);
    }
    PSRead(amrex::Vector<amrex::Geometry>& _geom, Unit a_unit,
           IO::ParmParse& pp, std::string name) : IC(_geom), field_unit(a_unit)
    {
        pp.queryclass(name, *this);
    }
    void Define() {

    };

    void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar)
    {
        const Set::Vector size = Set::Size(geom[lev]);
        const int ncomp = a_phi[lev]->nComp();
        const auto prob_lo = geom[lev].ProbLoArray();
        const auto cell_size = geom[lev].CellSizeArray();
        const amrex::GpuArray<int, AMREX_SPACEDIM> periodic = {
            AMREX_D_DECL(   geom[0].isPeriodic(0),
                            geom[0].isPeriodic(1),
                            geom[0].isPeriodic(2))};
        std::vector<Set::Scalar> coordinates(
            X.size() * AMREX_SPACEDIM);
        for (int n = 0; n < static_cast<int>(X.size()); ++n)
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                coordinates[n * AMREX_SPACEDIM + d] = X[n](d);
        amrex::Gpu::AsyncArray<Set::Scalar> X_device(
            coordinates.data(), coordinates.size());
        amrex::Gpu::AsyncArray<Set::Scalar> R_device(R.data(), R.size());
        const Set::Scalar* X = X_device.data();
        const Set::Scalar* R = R_device.data();
        const int sphere_count = static_cast<int>(this->X.size());
        const Set::Scalar eps = this->eps;
        const bool invert = this->invert;

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType type = a_phi[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x =
                    Set::Position(i, j, k, prob_lo, cell_size, type);
                Set::Scalar min_grain_id = 0;

                for (int n = 0; n < sphere_count; n++)
                {
                    Set::Scalar d;
                    Set::Scalar d1 = std::numeric_limits<Set::Scalar>::infinity();

                    Set::Vector center;
                    for (   int direction = 0;
                            direction < AMREX_SPACEDIM; ++direction)
                        center(direction) =
                            X[n * AMREX_SPACEDIM + direction];
                    d = (x - center).lpNorm<2>();

                    if (periodic[0])
                    {
                        d1 = std::min((x - center + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                            (x - center - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                    }
#if AMREX_SPACEDIM > 1
                    if (periodic[1])
                    {
                        d1 = std::min((x - center + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                            (x - center - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                    }
#endif
#if AMREX_SPACEDIM > 2
                    if (periodic[2])
                    {
                        d1 = std::min((x - center + size(2) * Set::Vector::Unit(2)).lpNorm<2>(),
                            (x - center - size(2) * Set::Vector::Unit(2)).lpNorm<2>());
                    }
#endif

                    d = std::min(d, d1);

                    if (d <= (R[n] + eps))

                    {
                        Set::Scalar m = 0.5 * (1 + erf((-d + R[n]) / (eps)));
                        min_grain_id = min_grain_id + m * (1. - min_grain_id);
                    }
                }

                Set::Scalar mask = 1.0;
                if (use_plane)
                {
                    const Set::Scalar distance =
                        (x - plane_point).dot(plane_normal);
                    mask = plane_eps > 0.0 ?
                        0.5 * (1.0 + std::tanh(
                            2.0 * distance / plane_eps)) :
                        static_cast<Set::Scalar>(distance >= 0.0);
                }

                phi(i, j, k, 0) =
                    field_value * mask *
                    (invert ? 1.0 - min_grain_id : min_grain_id);
                if (ncomp > 1)
                    phi(i, j, k, 1) =
                        field_value * mask * (1.0 - min_grain_id);
            });
        }
    }

    static void Parse(PSRead& value, IO::ParmParse& pp)
    {
        std::string filename;
        int verbose = 0;
        // Diffuseness of the sphere boundary
        pp.query_default("eps", value.eps, "0.0", Unit::Length());

        pp.forbid("filename", "use file.name instead");
        // Location of .xyzr file
        pp.query_file("file.name", filename); 

        // Verbosity (used in parser only)
        pp.query_default("verbose", verbose, 0); 

        // Unitless coordinate multiplier
        pp.query_default("mult",value.mult, "1.0", Unit::Less()); 

        // Whether to invert (1 outside instead of inside spheres)
        pp.query_default("invert",value.invert,false); 

        // Value represented by the packed-sphere indicator
        pp.query_default("value", value.field_value, "1.0", value.field_unit);

        // X-offset 
        pp.queryarr_default("x0",value.x0,"0 0 0", Unit::Length());

        if (pp.contains("x0"))
            pp_queryarr("x0",value.x0); // Coordinate offset

        // Optionally restrict the packing to one side of a diffuse plane.
        value.use_plane = pp.contains("plane.normal");
        if (value.use_plane)
        {
            pp.queryarr("plane.normal", value.plane_normal);
            pp.queryarr_default(
                "plane.point", value.plane_point, "0 0 0", Unit::Length());
            pp.query_default(
                "plane.eps", value.plane_eps, "0.0", Unit::Length());
            const Set::Scalar norm = value.plane_normal.norm();
            Util::AssertException(
                INFO, TEST(norm > 0.0), "plane.normal must be nonzero");
            Util::AssertException(
                INFO, TEST(value.plane_eps >= 0.0),
                "plane.eps must be nonnegative");
            value.plane_normal /= norm;
        }

        const bool use_slice = pp.contains("slice.z");
        Set::Scalar slice_z = 0.0;
        if (use_slice)
            pp.query_required("slice.z", slice_z, Unit::Length());
        
        Unit unit = Unit::Length();
        // Units of length in the file
        pp.queryunit("file.unit",unit);

        if (IO::ParmParse::InTraversalMode()) return;

        Util::AssertException(  INFO,TEST(unit.isType(Unit::Length())), 
                                "Unit must be of type length but got unit ",unit.normalized_unitstring());
        
        // Only the I/O processor touches the packing file.  This avoids
        // stressing a shared filesystem and, more importantly, guarantees
        // that every rank initializes from the exact same sphere list.
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            std::ifstream datafile(filename);
            std::string line;
            if (!datafile.is_open())
                Util::Abort(INFO, "Unable to open file ", filename);

            value.X.clear();
            value.R.clear();
            while (getline(datafile, line))
            {
                std::istringstream in(line);

                std::string strx, stry, strz, strR;
                in >> strx >> stry >> strz;
#if AMREX_SPACEDIM == 2
                // A 2-D packing conventionally stores x, y, and radius.
                // Accept that three-column form while retaining support for
                // x, y, z, radius files (the z coordinate is then ignored).
                if (!(in >> strR))
                {
                    strR = strz;
                    strz = "0";
                }
#else
                in >> strR;
#endif

                Set::Scalar x = (std::stod(strx) * unit).normalized_value();
                Set::Scalar y = (std::stod(stry) * unit).normalized_value();
                Set::Scalar z = (std::stod(strz) * unit).normalized_value();
                Set::Scalar r = (std::stod(strR) * unit).normalized_value();

#if AMREX_SPACEDIM == 2
                if (use_slice)
                {
                    const Set::Scalar dz = slice_z - z;
                    if (Util::Abs(dz) >= r) continue;
                    r = std::sqrt(r * r - dz * dz);
                }
#endif
                Set::Vector X(AMREX_D_DECL(x,y,z));
                X = value.x0 + value.mult*X;

                value.X.push_back(X);
                value.R.push_back(r);
                if (verbose > 0)
                    Util::Message(INFO, "x=", value.X.back().transpose(), " r=", value.R.back());
            }
            datafile.close();
        }

        int sphere_count = static_cast<int>(value.X.size());
        amrex::ParallelDescriptor::Bcast(
            &sphere_count, 1,
            amrex::ParallelDescriptor::IOProcessorNumber());

        std::vector<Set::Scalar> sphere_data(
            sphere_count * (AMREX_SPACEDIM + 1));
        if (amrex::ParallelDescriptor::IOProcessor())
        {
            for (int n = 0; n < sphere_count; ++n)
            {
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    sphere_data[n * (AMREX_SPACEDIM + 1) + d] = value.X[n](d);
                sphere_data[n * (AMREX_SPACEDIM + 1) + AMREX_SPACEDIM] =
                    value.R[n];
            }
        }

        if (!sphere_data.empty())
            amrex::ParallelDescriptor::Bcast(
                sphere_data.data(), static_cast<int>(sphere_data.size()),
                amrex::ParallelDescriptor::IOProcessorNumber());

        if (!amrex::ParallelDescriptor::IOProcessor())
        {
            value.X.resize(sphere_count);
            value.R.resize(sphere_count);
            for (int n = 0; n < sphere_count; ++n)
            {
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    value.X[n](d) =
                        sphere_data[n * (AMREX_SPACEDIM + 1) + d];
                value.R[n] =
                    sphere_data[n * (AMREX_SPACEDIM + 1) + AMREX_SPACEDIM];
            }
        }
    }

private:
    std::vector<Set::Vector> X;
    std::vector<Set::Scalar> R;
    Set::Scalar eps;
    Set::Scalar mult = 1.0;
    Set::Vector x0 = Set::Vector::Zero();
    bool invert = false;
    Unit field_unit = Unit::Less();
    Set::Scalar field_value = 1.0;
    bool use_plane = false;
    Set::Vector plane_point = Set::Vector::Zero();
    Set::Vector plane_normal = Set::Vector::Zero();
    Set::Scalar plane_eps = 0.0;
};
}
#endif
