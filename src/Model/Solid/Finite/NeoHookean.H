#ifndef MODEL_SOLID_FINITE_NEOHOOKEAN_H_
#define MODEL_SOLID_FINITE_NEOHOOKEAN_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "Unit/Unit.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
class NeoHookean : public Solid<Set::Sym::Major>
{
public:
    AMREX_GPU_HOST_DEVICE NeoHookean() {};
    AMREX_GPU_HOST_DEVICE NeoHookean(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    AMREX_GPU_HOST_DEVICE ~NeoHookean() = default;

    AMREX_GPU_HOST_DEVICE
    Set::Scalar W(const Set::Matrix& a_F) const
    {
#if AMREX_SPACEDIM==2
        const Set::Scalar J = a_F(0, 0) * a_F(1, 1) - a_F(0, 1) * a_F(1, 0);
        const Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        const Set::Scalar tr = a_F(0, 0) * a_F(0, 0)
                             + a_F(0, 1) * a_F(0, 1)
                             + a_F(1, 0) * a_F(1, 0)
                             + a_F(1, 1) * a_F(1, 1) + 1.0;
        Set::Scalar w = 0.0;
        w += 0.5 * mu * (tr / J23 - 3.);
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0);
        return w;
#elif AMREX_SPACEDIM==3
        Eigen::Matrix3d F = a_F;

        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        Set::Scalar w = 0.0;
        w += 0.5 * mu * ((F * F.transpose()).trace() / J23 - 3.);
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0);
        return w;
#endif
    }
    AMREX_GPU_HOST_DEVICE
    Set::Matrix DW(const Set::Matrix& a_F) const
    {
#if AMREX_SPACEDIM==2
        const Set::Scalar J = a_F(0, 0) * a_F(1, 1) - a_F(0, 1) * a_F(1, 0);
        const Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        const Set::Scalar tr = a_F(0, 0) * a_F(0, 0)
                             + a_F(0, 1) * a_F(0, 1)
                             + a_F(1, 0) * a_F(1, 0)
                             + a_F(1, 1) * a_F(1, 1) + 1.0;
        Set::Matrix FinvT;
        FinvT(0, 0) =  a_F(1, 1) / J;
        FinvT(0, 1) = -a_F(1, 0) / J;
        FinvT(1, 0) = -a_F(0, 1) / J;
        FinvT(1, 1) =  a_F(0, 0) / J;

        Set::Matrix r_dw;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
            {
                r_dw(i, j) = mu * (a_F(i, j) / J23 - tr * FinvT(i, j) / (3. * J23))
                           + kappa * (J - 1) * J * FinvT(i, j);
            }
        return r_dw;
#elif AMREX_SPACEDIM==3
        Eigen::Matrix3d F = a_F;

        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        // NOTE: the chained expression F.inverse().transpose() (a nested
        // Transpose<Inverse<Matrix3>> expression template) faults on the GPU
        // ("unspecified launch failure"). Materialize the inverse first.
        Eigen::Matrix3d Finv = F.inverse();
        Eigen::Matrix3d FinvT = Finv.transpose();

        Eigen::Matrix3d dw = Eigen::Matrix3d::Zero();

        dw += mu * (F / J23 - (F * F.transpose()).trace() * FinvT / (3. * J23));
        dw += kappa * (J - 1) * J * FinvT;

        return dw;
#endif
    }
    AMREX_GPU_HOST_DEVICE
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& a_F) const
    {
#if AMREX_SPACEDIM==2
        const Set::Scalar J = a_F(0, 0) * a_F(1, 1) - a_F(0, 1) * a_F(1, 0);
        const Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        const Set::Scalar tr = a_F(0, 0) * a_F(0, 0)
                             + a_F(0, 1) * a_F(0, 1)
                             + a_F(1, 0) * a_F(1, 0)
                             + a_F(1, 1) * a_F(1, 1) + 1.0;
        Set::Matrix FinvT;
        FinvT(0, 0) =  a_F(1, 1) / J;
        FinvT(0, 1) = -a_F(1, 0) / J;
        FinvT(1, 0) = -a_F(0, 1) / J;
        FinvT(1, 1) =  a_F(0, 0) / J;

        Set::Matrix4<2, Set::Sym::Major> r_ddw;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                for (int k = 0; k < 2; k++)
                    for (int l = 0; l < 2; l++)
                    {
                        Set::Scalar t1 = 0.0, t2 = 0.0;

                        if (i == k && j == l) t1 += 1.0;
                        t1 -= (2. / 3.) * a_F(i, j) * FinvT(k, l);
                        t1 -= (2. / 3.) * FinvT(i, j) * a_F(k, l);
                        t1 += (2. / 9.) * tr * FinvT(i, j) * FinvT(k, l);
                        t1 += (1. / 3.) * tr * FinvT(i, l) * FinvT(k, j);

                        t2 += (2. * J - 1.) * FinvT(i, j) * FinvT(k, l);
                        t2 += (1. - J) * FinvT(i, l) * FinvT(k, j);

                        r_ddw(i, j, k, l) = (mu / J23) * t1 + kappa * J * t2;
                    }
        return r_ddw;
#elif AMREX_SPACEDIM==3
        Eigen::Matrix3d F = a_F;

        Set::Matrix4<3, Set::Sym::Major> ddw;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        // NOTE: the chained expression F.inverse().transpose() (a nested
        // Transpose<Inverse<Matrix3>> expression template) faults on the GPU
        // ("unspecified launch failure"). Materialize the inverse first.
        Eigen::Matrix3d Finv = F.inverse();
        Eigen::Matrix3d FinvT = Finv.transpose();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        ddw(i, j, k, l) = 0.0;

                        Set::Scalar t1 = 0.0, t2 = 0.0;

                        if (i == k && j == l) t1 += 1.0;
                        t1 -= (2. / 3.) * F(i, j) * FinvT(k, l);
                        t1 -= (2. / 3.) * FinvT(i, j) * F(k, l);
                        t1 += (2. / 9.) * (F * F.transpose()).trace() * FinvT(i, j) * FinvT(k, l);
                        t1 += (1. / 3.) * (F * F.transpose()).trace() * FinvT(i, l) * FinvT(k, j);

                        t2 += (2. * J - 1.) * FinvT(i, j) * FinvT(k, l);
                        t2 += (1. - J) * FinvT(i, l) * FinvT(k, j);

                        ddw(i, j, k, l) = (mu / J23) * t1 + kappa * J * t2;
                    }
        return ddw;
#endif
    }
    AMREX_GPU_HOST_DEVICE
    void Print(std::ostream& out) const
    {
        AMREX_IF_ON_HOST(out << "mu = " << mu << " kappa = " << kappa;)
    }

public:
    Set::Scalar mu = NAN, kappa = NAN;
    static constexpr KinematicVariable kinvar = KinematicVariable::F;

public:
    static NeoHookean Zero()
    {
        NeoHookean ret;
        ret.mu = 0.0;
        ret.kappa = 0.0;
        return ret;
    }
    static NeoHookean Random()
    {
        NeoHookean ret;
        ret.mu = Util::Random();
        ret.kappa = Util::Random();
        return ret;
    }
    static void Parse(NeoHookean& value, IO::ParmParse& pp)
    {

        std::pair<std::string, Set::Scalar> moduli[2];

        pp.forbid("lame","Use 'lambda' instead for lame constant");
        pp.forbid("shear","Use 'mu' instead for shear modulus");
        pp.forbid("bulk","Use 'K' instead for bulk modulus");

        // Lame constant \(\lambda\),
        // shear modulus \(\mu\), Young's modulus \(E\), Poisson's ratio \(\nu\),
        // bulk modulus \(K\).
        // You can currently specify (mu and kappa), (lambda and mu), or (E and nu).
        // std::vector{Unit::Pressure, Unit::Pressure, Unit::Pressure, Unit::Less, Unit::Pressure};
        pp.query_exactly<2>({"lambda","mu","E","nu","kappa"}, moduli, std::vector<Unit>{Unit::Pressure(), Unit::Pressure(), Unit::Pressure(), Unit::Less(), Unit::Pressure()});


        if (moduli[0].first == "mu" && moduli[1].first == "kappa")
        {
            value.mu = moduli[0].second;
            value.kappa =  moduli[1].second;
        }
        else if (moduli[0].first == "lambda" && moduli[1].first == "mu")
        {
            Set::Scalar lambda = moduli[0].second;
            value.mu = moduli[1].second;
            value.kappa = lambda + (2.0 * value.mu) / 3.0;
        }
        else if (pp.contains("E") && pp.contains("nu")) {
            Set::Scalar E = moduli[0].second;
            Set::Scalar nu = moduli[1].second;
            value.kappa = E / (3.0 - 6.0 * nu);
            value.mu = E / (2.0 + 2.0 * nu);
        }
        else
        {
            Util::Exception(INFO,"Haven't implemented",moduli[0].first," and ",moduli[1].first," yet (sorry!)");
        }
    }

#define OP_CLASS NeoHookean
#define OP_VARS X(kappa) X(mu)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
