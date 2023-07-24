//
// Initialize a field using a mathematical expression.
// Expressions are imported as strings and are compiled real-time using the
// `AMReX Parser <https://amrex-codes.github.io/amrex/docs_html/Basics.html#parser>`_.
//
// Works for single or multiple-component fields.
// Use the :code:`regionN` (N=0,1,2, etc. up to number of components) to pass expression.
// For example:
//
// .. code-block:: bash
// 
//    ic.region0 = "sin(x*y*z)"
//    ic.region1 = "3.0*(x > 0.5 and y > 0.5)"
//
// for a two-component field. It is up to you to make sure your expressions are parsed
// correctly; otherwise you will get undefined behavior.
//


#ifndef IC_EXPRESSION_H_
#define IC_EXPRESSION_H_
#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "AMReX_Parser.H"

namespace IC
{
class Expression : public IC
{
private:
    enum CoordSys { Cartesian, Polar, Spherical };
    std::vector<amrex::Parser> parser;
    std::vector<amrex::ParserExecutor<4>> f;
    Expression::CoordSys coord = Expression::CoordSys::Cartesian;
public:
    Expression(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}
    Expression(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp.queryclass(name, *this);
    }
    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar a_time = 0.0)
    {
        Util::Assert(INFO, TEST(a_field[lev]->nComp() == (int)f.size()));
        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;// = mfi.tilebox();
            //bx.grow(a_field[lev]->nGrow());
            amrex::IndexType type = a_field[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            for (unsigned int n = 0; n < f.size(); n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                    if (coord == Expression::CoordSys::Cartesian)
                    {
#if AMREX_SPACEDIM == 1
                        field(i, j, k, n) = f[n](x(0), 0.0, 0.0, a_time);
#elif AMREX_SPACEDIM == 2
                        field(i, j, k, n) = f[n](x(0), x(1), 0.0, a_time);
#elif AMREX_SPACEDIM == 3
                        field(i, j, k, n) = f[n](x(0), x(1), x(2), a_time);
#endif
                    }
#if AMREX_SPACEDIM>1
                    else if (coord == Expression::CoordSys::Polar)
                    {
                        field(i, j, k, n) = f[n](sqrt(x(0)* x(0) + x(1) * x(1)), std::atan2(x(1), x(0)), x(2), a_time);
                    }
#endif
                });
            }
        }
        a_field[lev]->FillBoundary();
    };
    static void Parse(Expression& value, IO::ParmParse& pp)
    {
        for (int i = 0; true; i++)
        {
            std::string coordstr = "cartesian";
            // coordinate system to use: "cartesian" (for x,y,z,t) and 
            // "polar" (for r, theta, z, t)
            pp.query("coord", coordstr); 
            if (coordstr == "cartesian") value.coord = Expression::CoordSys::Cartesian;
            else if (coordstr == "polar") value.coord = Expression::CoordSys::Polar;
            else Util::Abort(INFO, "unsupported coordinates ", coordstr);

            std::string func = "0.0";
            std::string name = "region" + std::to_string(i);

            if (!pp.contains(name.data())) break;
            pp.query(name.data(), func);

            value.parser.push_back(amrex::Parser(func));

            if (value.coord == Expression::CoordSys::Cartesian)
            {
                value.parser.back().registerVariables({ "x","y","z","t" });
                value.f.push_back(value.parser.back().compile<4>());
            }
            else if (value.coord == Expression::CoordSys::Polar)
            {
                value.parser.back().registerVariables({ "r","theta","z","t" });
                value.f.push_back(value.parser.back().compile<4>());
            }
        }

    };
};
}

#endif
