#ifndef MODEL_SOLID_AFFINE_ISOTROPIC_H_
#define MODEL_SOLID_AFFINE_ISOTROPIC_H_

#include "AMReX.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Affine/Affine.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class Isotropic : public Affine<Set::Sym::Isotropic>
{
public:

    Isotropic() {};
    Isotropic(Affine<Set::Sym::Isotropic> base) : Affine<Set::Sym::Isotropic>(base) {};
    Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0 = Set::Matrix::Zero())
    {
        Define(a_mu, a_lambda, a_F0);
    };

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0)
    {
        F0 = a_F0;
        ddw = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic>(a_lambda, a_mu);
    }

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
    {
        ddw = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic>(a_lambda, a_mu);
    }

    Set::Scalar W(const Set::Matrix& F) const override
    {
        return 0.5 * ((F - F0).transpose() * (ddw * ((F - F0)))).trace();
    }
    Set::Matrix DW(const Set::Matrix& F) const override
    {
        return ddw * (F - F0);
    }
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Isotropic> DDW(const Set::Matrix& /*F*/) const override
    {
        return ddw;
    }
    virtual void Print(std::ostream& out) const override
    {
        out << "lambda=" << ddw.Lambda() << ", mu=" << ddw.Mu() << ", F0=\n" << F0;
    }
    virtual bool ContainsNan() override
    {
        if (F0.hasNaN()) return true;
        if (ddw.contains_nan()) return true;
        else return false;
    }
public:
    static Isotropic Zero()
    {
        Isotropic ret;
        Set::Scalar mu = 0.0;
        Set::Scalar lambda = 0.0;
        Set::Matrix F0 = Set::Matrix::Zero();
        ret.Define(mu, lambda, F0);
        return ret;
    }
    static Isotropic Random()
    {
        Isotropic ret;
        Set::Scalar mu = Util::Random();
        Set::Scalar lambda = Util::Random();
        Set::Matrix F0 = Set::Matrix::Random();
        ret.Define(mu, lambda, F0);
        return ret;
    }
    static void Parse(Isotropic& value, IO::ParmParse& pp)
    {
        Set::Scalar mu = NAN, lambda = NAN;
        Set::Matrix F0 = Set::Matrix::Zero();

        std::pair<std::string, Set::Scalar> moduli[2];

        pp.forbid("lame","Use 'lambda' instead for lame constant");
        pp.forbid("shear","Use 'mu' instead for shear modulus");
        pp.forbid("bulk","Use 'kappa' instead for bulk modulus");

        // Specify exactly two of: lame constant \(\lambda\),
        // shear modulus \(\mu\), Young's modulus \(E\), Poisson's ratio \(\nu\),
        // bulk modulus \(\kappa\).
        // \(\mu\) and \(\lambda\) are how the final values are stored.
        pp.query_exactly<2>({"lambda","mu","E","nu","kappa"}, moduli);

        if      (moduli[0].first == "lambda" && moduli[1].first == "mu")
        {
            lambda = moduli[0].second;
            mu = moduli[1].second;
        }
        else if (moduli[0].first == "lambda" && moduli[1].first == "E")
        {
            lambda = moduli[0].second;
            Set::Scalar E = moduli[1].second;
            mu = (E - 3.0 * lambda + sqrt(E * E + 9.0 * lambda * lambda + 2.0 * E * lambda)) / 4.0;
        }
        else if (moduli[0].first == "lambda" && moduli[1].first == "nu")
        {
            lambda = moduli[0].second;
            Set::Scalar nu = moduli[1].second;
            mu = lambda * (1.0 - 2.0 * nu) / 2.0 / nu;
        }
        else if (moduli[0].first == "mu" && moduli[1].first == "E")
        {
            mu = moduli[0].second;
            Set::Scalar E = moduli[1].second;
            lambda = mu * (E - 2.0 * mu) / (3.0 * mu - E);
        }
        else if (moduli[0].first == "mu" && moduli[1].first == "nu")
        {
            mu = moduli[0].second;
            Set::Scalar nu = moduli[1].second;
            lambda = 2.0 * mu * nu / (1.0 - 2.0 * nu);
        }
        else if (moduli[0].first == "mu" && moduli[1].first == "K")
        {
            mu = moduli[0].second;
            Set::Scalar K = moduli[1].second;
            lambda = K - (2.0 * mu / 3.0);
        }
        else if (moduli[0].first == "E" && moduli[1].first == "nu")
        {
            Set::Scalar E = moduli[0].second, nu = moduli[1].second;
            lambda = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu);
            mu = E / 2.0 / (1.0 + nu);
        }
        else
        {
            Util::Exception(INFO,"Haven't implemented ",moduli[0].first," and ",moduli[1].first," (sorry!)");
        }

        if (pp.contains("F0"))
        {
            pp_queryarr("F0", F0); // Eigendeformation gradient
        }
        value.Define(mu, lambda, F0);
    }
#define OP_CLASS Isotropic
#define OP_VARS  X(ddw) X(F0)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"
}
}
}

template<>
ALAMO_SINGLE_DEFINITION
int Set::Field<Model::Solid::Affine::Isotropic>::NComp() const
{
    return AMREX_SPACEDIM * AMREX_SPACEDIM + 2;
}
template<>
ALAMO_SINGLE_DEFINITION
std::string Set::Field<Model::Solid::Affine::Isotropic>::Name(int i) const
{
    if (i == 0) return name + "_lambda";
    if (i == 1) return name + "_mu";
#if AMREX_SPACEDIM==2
    if (i == 2) return name + "_F0_xx";
    if (i == 3) return name + "_F0_xy";
    if (i == 4) return name + "_F0_yx";
    if (i == 5) return name + "_F0_yy";
#elif AMREX_SPACEDIM==3
    if (i == 2) return name + "_F0_xx";
    if (i == 3) return name + "_F0_xy";
    if (i == 4) return name + "_F0_xy";
    if (i == 5) return name + "_F0_yx";
    if (i == 6) return name + "_F0_yy";
    if (i == 7) return name + "_F0_yz";
    if (i == 8) return name + "_F0_zy";
    if (i == 9) return name + "_F0_zx";
    if (i == 10) return name + "_F0_zz";
#endif
    return name;
}
template<>
ALAMO_SINGLE_DEFINITION
void Set::Field<Model::Solid::Affine::Isotropic>::Copy(int a_lev, amrex::MultiFab& a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Affine::Isotropic> const& src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const& dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM * AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i, j, k, a_dstcomp + 0) = src(i, j, k).ddw.Lambda();
                    dst(i, j, k, a_dstcomp + 1) = src(i, j, k).ddw.Mu();
#if AMREX_SPACEDIM==2
                    dst(i, j, k, a_dstcomp + 2) = src(i, j, k).F0(0, 0);
                    dst(i, j, k, a_dstcomp + 3) = src(i, j, k).F0(0, 1);
                    dst(i, j, k, a_dstcomp + 4) = src(i, j, k).F0(1, 0);
                    dst(i, j, k, a_dstcomp + 5) = src(i, j, k).F0(1, 1);
#elif AMREX_SPACEDIM==3
                    dst(i, j, k, a_dstcomp + 2) = src(i, j, k).F0(0, 0);
                    dst(i, j, k, a_dstcomp + 3) = src(i, j, k).F0(0, 1);
                    dst(i, j, k, a_dstcomp + 4) = src(i, j, k).F0(0, 2);
                    dst(i, j, k, a_dstcomp + 5) = src(i, j, k).F0(1, 0);
                    dst(i, j, k, a_dstcomp + 6) = src(i, j, k).F0(1, 1);
                    dst(i, j, k, a_dstcomp + 7) = src(i, j, k).F0(1, 2);
                    dst(i, j, k, a_dstcomp + 8) = src(i, j, k).F0(2, 0);
                    dst(i, j, k, a_dstcomp + 9) = src(i, j, k).F0(2, 1);
                    dst(i, j, k, a_dstcomp + 10) = src(i, j, k).F0(2, 2);
#endif
                });
            }
        }
    }
}




#endif
