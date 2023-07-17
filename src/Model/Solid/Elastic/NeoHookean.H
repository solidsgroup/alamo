#ifndef MODEL_SOLID_ELASTIC_NEOHOOKEAN_H_
#define MODEL_SOLID_ELASTIC_NEOHOOKEAN_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class NeoHookean : public Solid<Set::Sym::Major>
{
public:
    NeoHookean() {};
    NeoHookean(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~NeoHookean() {};

    Set::Scalar W(const Set::Matrix& F) const override
    {
        Set::Matrix Ft = F * F0.inverse();
        Set::Scalar J = Ft.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        Set::Scalar w = 0.0;
        w += 0.5 * mu * ((Ft * Ft.transpose()).trace() / J23 - 3.);
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0);
        return w;
    }
    Set::Matrix DW(const Set::Matrix& F) const override
    {
        Set::Matrix Ft = F * F0.inverse();
        Set::Scalar J = Ft.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        Set::Matrix FinvT = Ft.inverse().transpose();

        Set::Matrix dw = Set::Matrix::Zero();

        dw += mu * (Ft / J23 - (Ft * Ft.transpose()).trace() * FinvT / (3. * J23));
        dw += kappa * (J - 1) * J * FinvT;

        return dw;
    }
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& F) const override
    {
        Set::Matrix Ft = F * F0.inverse();
        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> ddw;
        Set::Scalar J = Ft.determinant();
        Set::Scalar J23 = std::pow(fabs(J), 2. / 3.);
        Set::Matrix FinvT = Ft.inverse().transpose();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        ddw(i, j, k, l) = 0.0;

                        Set::Scalar t1 = 0.0, t2 = 0.0;

                        if (i == k && j == l) t1 += 1.0;
                        t1 -= (2. / 3.) * Ft(i, j) * FinvT(k, l);
                        t1 -= (2. / 3.) * FinvT(i, j) * F(k, l);
                        t1 += (2. / 9.) * (Ft * Ft.transpose()).trace() * FinvT(i, j) * FinvT(k, l);
                        t1 += (1. / 3.) * (Ft * Ft.transpose()).trace() * FinvT(i, l) * FinvT(k, j);

                        t2 += (2. * J - 1.) * FinvT(i, j) * FinvT(k, l);
                        t2 += (1. - J) * FinvT(i, l) * FinvT(k, j);

                        ddw(i, j, k, l) = (mu / J23) * t1 + kappa * J * t2;
                    }
        return ddw;
    }
    virtual void Print(std::ostream& out) const override
    {
        out << "mu = " << mu << " kappa = " << kappa;
    }

public:
    Set::Scalar mu = NAN, kappa = NAN;
    Set::Matrix F0 = Set::Matrix::Zero();
    static constexpr KinematicVariable kinvar = KinematicVariable::F;

public:
    static NeoHookean Zero()
    {
        NeoHookean ret;
        ret.mu = 0.0;
        ret.kappa = 0.0;
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }
    static NeoHookean Random()
    {
        NeoHookean ret;
        ret.mu = Util::Random();
        ret.kappa = Util::Random();
        ret.F0 = Set::Matrix::Random();
        return ret;
    }
    static void Parse(NeoHookean& value, IO::ParmParse& pp)
    {
        //Set::Scalar mu = NAN, kappa = NAN;
        if (pp.contains("shear") && pp.contains("kappa")) {
            pp.query("shear", value.mu);
            pp.query("kappa", value.kappa);
        }
        else if (pp.contains("lame") && pp.contains("shear")) {
            Set::Scalar lame;
            pp.query("shear", value.mu);
            pp.query("lame", lame);
            value.kappa = lame + (2.0 * value.mu) / 3.0;
        }
        else if (pp.contains("E") && pp.contains("nu")) {
            Set::Scalar E, nu;
            pp.query("E", E);
            pp.query("nu", nu);
            value.kappa = E / (3.0 - 6.0 * nu);
            value.mu = E / (2.0 + 2.0 * nu);
        }
        if (pp.contains("F0")) {
            pp.queryarr("F0", value.F0);
        }
    }

#define OP_CLASS NeoHookean
#define OP_VARS X(kappa) X(mu) X(F0)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
