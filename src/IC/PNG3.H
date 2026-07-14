//
// Initialize a 3-component field using RGB channels from a PNG image. (2D only)
//

#ifndef IC_PNG3_H
#define IC_PNG3_H
#include <cmath>

#ifndef ALAMO_NOPNG
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <png.h>
#endif

#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"
#include "IO/ParmParse.H"
#include "Util/PNG.H"

namespace IC
{
class PNG3 : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "png3";

    virtual ~PNG3() {}

    PNG3(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PNG3(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : PNG3(_geom)
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
                if (type == amrex::IndexType::TheNodeType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i)) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j)) * geom[lev].CellSize()[1];
                }
                else if (type == amrex::IndexType::TheCellType())
                {
                    x(0) = domlo(0) + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
                    x(1) = domlo(1) + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
                }

                std::array<Set::Scalar, 4> val = png(x);

                field(i, j, k, 0) = val[0];  // R
                field(i, j, k, 1) = val[1];  // G
                field(i, j, k, 2) = val[2];  // B
            });
        }
        a_field[lev]->FillBoundary();
#else
        Util::IgnoreUnused(lev, a_field);
        Util::Abort(INFO, "PNG is disabled");
#endif
    };

private:
#ifndef ALAMO_NOPNG
    Util::PNG png;
#endif

public:
    static void Parse(PNG3& value, IO::ParmParse& pp)
    {
#ifndef ALAMO_NOPNG
        pp.queryclass<Util::PNG>(value.png);
#else
        Util::IgnoreUnused(value, pp);
        Util::Abort(INFO, "PNG is disabled");
#endif
    }
};
}
#endif