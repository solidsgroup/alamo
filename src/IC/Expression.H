#ifndef IC_EXPRESSION_H_
#define IC_EXPRESSION_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

#include "libmesh/fparser_ad.hh"

/// \class Expression
/// \brief Use JIT math expression parsing to read arbitrary domain description
namespace IC
{
    class Expression : public IC
    {
    public:
        Expression(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}

        void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
        {
            for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(a_field[lev]->nGrow());
                amrex::IndexType type = a_field[lev]->ixType();

                amrex::Array4<Set::Scalar> const &field = a_field[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                    
                    //Set::Scalar *loc = x.data();
                    double loc[3];
                    loc[0] = x(0); loc[1] = x(1); loc[2] = x(2);

                    field(i,j,k,0) = fparser[0].Eval(loc);
                    field(i,j,k,1) = fparser[1].Eval(loc);
                    field(i,j,k,2) = fparser[2].Eval(loc);

//                    Set::Scalar min_value = 100.0;
//                    for (int i = 0; i < center.size(); i++)
//                    {
//                        Set::Scalar value = 0.;
//                        Set::Scalar norm = (A[i] * (x - center[i])).lpNorm<2>();
//                        value = 0.5 + 0.5 * std::erf(((x - center[i]).transpose() * A[i] * (x - center[i]) - 1.0) / eps[i] / norm);
//                        min_value = value < min_value ? value : min_value;
//                    }
//                    field(i, j, k) = min_value;
//                    if (field(i, j, k) < 0.)
//                        field(i, j, k) = 0.;
//                    if (field(i, j, k) > 1.)
//                        field(i, j, k) = 1.;
                });
            }
            a_field[lev]->FillBoundary();
        };

    private:
    //
    // TODO: this part still needs to be updated. Cannot yet handle arbitrary number of grains...
    //
        std::array<FunctionParserAD,3> fparser;
        //	int dim = AMREX_SPACEDIM;
        //	amrex::Vector<Set::Vector> center;
        //	amrex::Vector<Set::Scalar> eps;
        //	amrex::Vector<Set::Matrix> A;

    public:
        static void Parse(Expression &value, IO::ParmParse &pp)
        {
            std::array<char *,3> names = {"region1","region2","region3"};
            for (int i = 0; i < value.fparser.size(); i++)
            {
                std::string func = "0.0";
                pp.query(names[i],func);

                value.fparser[i].AddVariable("x");
                value.fparser[i].AddVariable("y");
                value.fparser[i].AddVariable("z");
                if (value.fparser[i].Parse(func.c_str(), "x,y,z") != -1)
                {
                    std::cout << "Parsing failed for " << func << '\n';
                    exit(1);
                }
                value.fparser[i].Optimize();
                if (!value.fparser[i].JITCompile())
                {
                    std::cout << "JIT compile failed for " << func << '\n';
                    exit(1);
                }
            }
        }
    };
} // namespace IC
#endif
