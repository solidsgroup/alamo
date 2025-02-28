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

                field(i,j,k,0) = val[channel];


                if (field.nComp() > 1) field(i, j, k, 1) = 1.0 - field(i, j, k, 0);

            });
        }
        a_field[lev]->FillBoundary();
#else
        Util::Abort(INFO,"PNG is disabled");
#endif
        
    };

private:
#ifndef ALAMO_NOPNG
    Util::PNG png;
    Channel channel = Channel::G;
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

        
        pp.queryclass<Util::PNG>(value.png);

#else
        Util::Abort(INFO,"PNG is disabled");
#endif
    }
};
}
#endif
