#ifndef MODEL_SOLID_AFFINE_CUBIC_H_
#define MODEL_SOLID_AFFINE_CUBIC_H_

#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/Cubic.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class Cubic : public Linear::Cubic
{
public:

	Cubic() {};
    Cubic(Solid<Set::Sym::MajorMinor> base) : Linear::Cubic(base) {};
	virtual ~Cubic() {};

    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2,Set::Matrix a_F0=Set::Matrix::Zero())
    {
        Linear::Cubic::Define(C11,C12,C44,phi1,Phi,phi2);
        F0 = a_F0;
    }
    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R, Set::Matrix a_F0=Set::Matrix::Zero())
    {
        Linear::Cubic::Define(C11,C12,C44,R);
        F0 = a_F0;
    }
    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return Linear::Cubic::W(gradu - F0);
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return Linear::Cubic::DW(gradu - F0);
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix & gradu) const override
    {
        return Linear::Cubic::DDW(gradu - F0);
    }
    virtual void Print(std::ostream &out) const //override 
    {
        out << ddw;
    }

public:
    Set::Matrix F0 = Set::Matrix::Zero();

    static Cubic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static Cubic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        Cubic ret = Linear::Cubic::Random(C11,C12,C44);
        ret.F0 = Set::Matrix::Random();
        return ret;
    }
    static void Parse(Cubic & value, IO::ParmParse & pp)
    {
        Linear::Cubic::Parse(value,pp);
        if (pp.contains("F0")) pp.queryarr("F0",value.F0);
    }

    AMREX_FORCE_INLINE
    static Cubic Combine(const std::vector<Cubic> &models, const std::vector<Set::Scalar> &eta)
    {
        Cubic ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        ret.F0 = Set::Matrix::Zero();
        Set::Scalar etasum = 0.;
        for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        for (unsigned int n = 0 ; n < models.size(); n++)
        {
            ret.ddw += models[n].ddw * (eta[n] / etasum);
            ret.F0  += models[n].F0  * (eta[n] / etasum);
        }
        //if (eta[0] > eta[1])
        //{
        //    ret.ddw = models[0].ddw;
        //    ret.F0  = models[0].F0;
        //}
        //else
        //{
        //    ret.ddw = models[1].ddw;
        //    ret.F0  = models[1].F0;
        //}
        return ret;
    }
};


}
}
}

#endif

