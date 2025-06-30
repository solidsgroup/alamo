#ifndef IC_BETA_DISTRIBUTION_H_
#define IC_BETA_DISTRIBUTION_H_

#include "IC/Random.H"
#include "IO/ParmParse.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "Util/Distribution.H"
#include "Util/Util.H"

#include "AMReX_Box.H"
#include "AMReX_Geometry.H"
#include "AMReX_GpuLaunchFunctsC.H"
#include "AMReX_GpuQualifiers.H"

#include "AMReX_IndexType.H"
#include "AMReX_MFIter.H"
#include "AMReX_Vector.H"

namespace IC
{
/// Set each point to a random value with a specified beta distribution
class BetaDistribution : public Random
{
public:
    static constexpr const char *name = "betadistribution";

    BetaDistribution(amrex::Vector<amrex::Geometry> &_geom, Set::Scalar a_mult)
        : Random(_geom, a_mult) {}
    BetaDistribution(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name)
        : Random(_geom, 1.0)
    {
        pp.queryclass(name, *this);
    }

    void
    Add(const int &lev, Set::Field<Set::Scalar> &a_field, Set::Scalar)
    {
        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType type = a_field[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())
                bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType())
                bx = mfi.grownnodaltilebox();
            else
                Util::Abort(INFO, "Unknown index type");

            Set::Patch<Set::Scalar> field = a_field.Patch(lev, mfi);

            // Capture parameters for GPU kernel
            Set::Scalar local_alpha = alpha;
            Set::Scalar local_beta = beta;
            Set::Scalar local_offset = offset;
            Set::Scalar local_mult = mult;

            for (int n = 0; n < a_field[lev]->nComp(); n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Scalar beta_value = Util::Distribution::Beta(local_alpha, local_beta);
                    field(i, j, k, n) += local_offset + local_mult * beta_value;
                });
            }
        }
    }

    using IC::Add;

    static void
    Parse(BetaDistribution &value, IO::ParmParse &pp)
    {
        Random::Parse(value, pp);

        // alpha parameter for a beta distribution
        pp.query_required("alpha", value.alpha);
        // beta parameter for a beta distribution
        pp.query_required("beta", value.beta);
        if (value.alpha <= 0.0 || value.beta <= 0.0)
            Util::Abort(INFO, "Beta distribution parameters alpha and beta must be positive");
    }

private:
    Set::Scalar alpha = NAN;
    Set::Scalar beta = NAN;
};
}
#endif
