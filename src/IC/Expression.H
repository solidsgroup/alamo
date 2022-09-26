#ifndef IC_EXPRESSION_H_
#define IC_EXPRESSION_H_
#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

#ifndef ALAMO_JIT
namespace IC 
{
    class Expression: public IC
    {
        public:
        Expression(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {Util::Abort(INFO,"To use IC::Expression, you must compile with the JIT option!");}
        Expression(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : IC(_geom)
        {pp.queryclass(name,*this);}
        void Add(const int &, Set::Field<Set::Scalar> &) {};
        static void Parse(Expression &, IO::ParmParse &) {};
    };
}
#else

#include "libmesh/fparser_ad.hh"
/// \class Expression
/// \brief Use JIT math expression parsing to read arbitrary domain description
namespace IC
{
    class Expression : public IC
    {
    public:
        enum CoordSys {Cartesian, Polar, Spherical};
        Expression(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
        Expression(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : IC(_geom) 
        {pp.queryclass(name,*this);}

        void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
        {
            Util::Assert(INFO,TEST(a_field[lev]->nComp() == (int)fparser.size()));
            for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(a_field[lev]->nGrow());
                amrex::IndexType type = a_field[lev]->ixType();

                amrex::Array4<Set::Scalar> const &field = a_field[lev]->array(mfi);
                for (unsigned int n = 0; n < fparser.size(); n++)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                        double loc[AMREX_SPACEDIM];
                        if (coord == Expression::CoordSys::Cartesian)
                        {
                            AMREX_D_TERM(loc[0] = x(0);, loc[1] = x(1);, loc[2] = x(2););
                            field(i,j,k,n) = fparser[n].Eval(loc);
                        }
                        #if AMREX_SPACEDIM>1
                        else if (coord == Expression::CoordSys::Polar)
                        {
                            AMREX_D_TERM(loc[0] = sqrt(x(0)*x(0) + x(1)*x(1));, loc[1] = std::atan2(x(1),x(0));, loc[2] = x(2));
                            field(i,j,k,n) = fparser[n].Eval(loc);
                        }
                        #endif
                    });
                }
            }
            a_field[lev]->FillBoundary();
        };

    private:
        std::vector<FunctionParserAD> fparser;
        Expression::CoordSys coord = Expression::CoordSys::Cartesian;

    public:
        // This is an extremely general IC that uses JIT fparser to compile an
        // IC based on a mathematical expression at runtime. Eventually we will
        // probably replace all ICs with this one. 
        // 
        // This can be used with an arbitrary number of components, named `region1`,
        // `region2`, etc. You do not need to specify a number ahead of time.
        // Each region should be a string that represents a function in terms of 
        // x, y, z, and t. 
        // It can be a boolean expression (returning 1 or 0) or it can return a 
        // value. 
        static void Parse(Expression &value, IO::ParmParse &pp)
        {
            for (int i = 0; true; i++)
            {
                std::string coordstr = "cartesian";
                pp.query("coord",coordstr);
                if (coordstr == "cartesian") value.coord = Expression::CoordSys::Cartesian;
                else if (coordstr == "polar") value.coord= Expression::CoordSys::Polar;
                else Util::Abort(INFO,"unsupported coordinates ",coordstr);

                std::string func = "0.0";
                std::string name = "region" + std::to_string(i);
                
                if (!pp.contains(name.data())) break;
                pp.query(name.data(),func);

                FunctionParserAD fparser;

                if (value.coord == Expression::CoordSys::Cartesian)
                {
                    fparser.AddVariable("x");
                    fparser.AddVariable("y");
                    fparser.AddVariable("z");
                    if (fparser.Parse(func.c_str(), "x,y,z") != -1)
                    {
                        std::cout << "Parsing failed for " << func << '\n';
                        exit(1);
                    }
                }
                else if (value.coord == Expression::CoordSys::Polar)
                {
                    fparser.AddVariable("r");
                    fparser.AddVariable("theta");
                    fparser.AddVariable("z");
                    if (fparser.Parse(func.c_str(), "r,theta,z") != -1)
                    {
                        std::cout << "Parsing failed for " << func << '\n';
                        exit(1);
                    }
                }

                fparser.Optimize();
                if (!fparser.JITCompile())
                {
                    std::cout << "JIT compile failed for " << func << '\n';
                    exit(1);
                }
                value.fparser.push_back(fparser);
            }
        }
    };
} // namespace IC
#endif
#endif
