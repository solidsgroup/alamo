// 
// Initialize a field using a PNG image. (2D only)
//

#ifndef IC_PNG_H
#define IC_PNG_H
#include <cmath>

#ifndef ALAMO_NOPNG
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <png.h>
#endif

#include "IC/IC.H"
#include "Util/Util.H"
#include "Util/BMP.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"
#include "Util/PNG.H"

namespace IC
{
class PNG : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "png";

    //enum Type {XYZ, XY, YZ, XZ};
    enum Fit { Stretch, FitWidth, FitHeight, Coord };
    enum Channel { R = 0, G = 1, B = 2, A = 3 };

    virtual ~PNG() {}

    PNG(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PNG(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : PNG(_geom)
    {
        pp_queryclass(name, *this);
    }

    PNG(amrex::Vector<amrex::Geometry>& _geom, Unit a_unit, IO::ParmParse& pp, std::string name) : PNG(_geom)
    {
        unit = a_unit;
        pp_queryclass(name, *this);
    }

    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar)
    {

#ifndef ALAMO_NOPNG
        Set::Vector DX(geom[lev].CellSize());
        amrex::Box domain = geom[lev].Domain();

        amrex::IndexType type = a_field[lev]->ixType();
        domain.convert(type);

        Set::Vector domlo(AMREX_D_DECL(geom[lev].ProbLo()[0], geom[lev].ProbLo()[1], 0.0));
        Set::Vector domhi(AMREX_D_DECL(geom[lev].ProbHi()[0], geom[lev].ProbHi()[1], 0.0));

        png.setDomain(domlo, domhi);
        const Channel selected_channel = channel;
        const Set::Scalar output_value = field_value;
        const bool invert_image = invert;
        const bool threshold_image = use_threshold;
        const Set::Scalar image_threshold = threshold;
        const bool mask_with_plane = use_plane;
        const Set::Vector mask_point = plane_point;
        const Set::Vector mask_normal = plane_normal;
        const Set::Scalar mask_eps = plane_eps;

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            if (type == amrex::IndexType::TheCellType()) bx = mfi.growntilebox();
            bx = bx & domain;

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x = Set::Vector::Zero();
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i)) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)) * geom[lev].CellSize()[1];
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i)+0.5) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)+0.5) * geom[lev].CellSize()[1];
                }


                std::array<Set::Scalar, 4> val = png(x);
                Set::Scalar indicator = val[selected_channel];
                if (threshold_image)
                    indicator = indicator >= image_threshold ? 1.0 : 0.0;
                if (invert_image) indicator = 1.0 - indicator;

                Set::Scalar mask = 1.0;
                if (mask_with_plane)
                {
                    const Set::Scalar distance =
                        (x - mask_point).dot(mask_normal);
                    mask = mask_eps > 0.0 ?
                        0.5 * (1.0 + std::tanh(
                            2.0 * distance / mask_eps)) :
                        static_cast<Set::Scalar>(distance >= 0.0);
                }

                field(i,j,k,0) = output_value * mask * indicator;
                if (field.nComp() > 1)
                    field(i,j,k,1) =
                        output_value * mask * (1.0 - indicator);

            });
        }
        a_field[lev]->FillBoundary();
#else
        Util::IgnoreUnused(lev,a_field);
        Util::Abort(INFO,"PNG is disabled");
#endif
        
    };

private:
#ifndef ALAMO_NOPNG
    Util::PNG png;
    Channel channel = Channel::G;
    Unit unit = Unit::Less();
    Set::Scalar field_value = 1.0;
    bool invert = false;
    bool use_threshold = false;
    Set::Scalar threshold = 0.5;
    bool use_plane = false;
    Set::Vector plane_point = Set::Vector::Zero();
    Set::Vector plane_normal = Set::Vector::Zero();
    Set::Scalar plane_eps = 0.0;
#endif

public:
    static void Parse(PNG& value, IO::ParmParse& pp)
    {
#ifndef ALAMO_NOPNG

        std::string channel = "g";
        // Color channel to use (options: r, R, g, G, b, B, a, A)        
        pp_query_validate("channel", channel, {"r","g","b","a","R","G","B","A"}); 
        if (channel == "r" || channel == "R") value.channel = Channel::R;
        else if (channel == "g" || channel == "G") value.channel = Channel::G;
        else if (channel == "b" || channel == "B") value.channel = Channel::B;
        else if (channel == "a" || channel == "A") value.channel = Channel::A;
        else Util::Abort(INFO, "Invalid value for bmp channel - should be r/g/b/a but received '", channel, "'");

        pp.query_default(
            "value", value.field_value, "1.0", value.unit);
        pp.query_default("invert", value.invert, false);
        value.use_threshold = pp.contains("threshold");
        if (value.use_threshold)
        {
            pp.query_required(
                "threshold", value.threshold, Unit::Less());
            if (value.threshold < 0.0 || value.threshold > 1.0)
                Util::Exception(
                    INFO, "PNG threshold must lie in [0,1]");
        }

        value.use_plane = pp.contains("plane.normal");
        if (value.use_plane)
        {
            pp.queryarr("plane.normal", value.plane_normal);
            pp.queryarr_default(
                "plane.point", value.plane_point,
                "0 0 0", Unit::Length());
            pp.query_default(
                "plane.eps", value.plane_eps,
                "0.0", Unit::Length());
            const Set::Scalar norm = value.plane_normal.norm();
            Util::AssertException(
                INFO, TEST(norm > 0.0),
                "PNG plane.normal must be nonzero");
            Util::AssertException(
                INFO, TEST(value.plane_eps >= 0.0),
                "PNG plane.eps must be nonnegative");
            value.plane_normal /= norm;
        }

        pp.queryclass<Util::PNG>(value.png);

#else
        if (IO::ParmParse::InTraversalMode()) return;
        Util::Abort(INFO,"PNG is disabled");
#endif
    }
};
}
#endif
