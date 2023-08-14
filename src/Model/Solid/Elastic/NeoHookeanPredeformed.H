#ifndef MODEL_SOLID_ELASTIC_NEOHOOKEANPREDEFORMED_H_
#define MODEL_SOLID_ELASTIC_NEOHOOKEANPREDEFORMED_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Elastic/NeoHookean.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class NeoHookeanPredeformed : public NeoHookean
{
public:
    NeoHookeanPredeformed() {};
    NeoHookeanPredeformed(NeoHookean base) : NeoHookean(base) {};
    virtual ~NeoHookeanPredeformed() {};

    Set::Scalar W(const Set::Matrix& F) const override
    {
        return NeoHookean::W(F*F0.inverse());
    }
    Set::Matrix DW(const Set::Matrix& F) const override
    {
        return NeoHookean::DW(F*F0.inverse());
    }
    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix& F) const override
    {
        return NeoHookean::DDW(F*F0.inverse());
    }
    virtual void Print(std::ostream& out) const override
    {
        out << "mu = " << mu << " kappa = " << kappa << " F0 = " << F0;
    }

public:
    Set::Matrix F0 = Set::Matrix::Identity();

public:
    static NeoHookeanPredeformed Zero()
    {
        NeoHookeanPredeformed ret = NeoHookean::Zero();
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }
    static NeoHookean Random()
    {
        NeoHookeanPredeformed ret = NeoHookean::Random();
        ret.F0 = Set::Matrix::Random();
        return ret;
    }
    static void Parse(NeoHookeanPredeformed& value, IO::ParmParse& pp)
    {
        NeoHookean::Parse(value,pp);
        if (pp.contains("F0")) {
            pp.queryarr("F0", value.F0);
        }
    }

#define OP_CLASS NeoHookeanPredeformed
#define OP_VARS X(kappa) X(mu) X(F0)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
