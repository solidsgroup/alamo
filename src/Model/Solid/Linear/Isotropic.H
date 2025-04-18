//
// This model implements an isotropic linear elastic material.
// See `this link <https://en.wikipedia.org/wiki/Linear_elasticity#(An)isotropic_(in)homogeneous_media>`_
// for more information about the theory.
//
// Free energy for a linear material is defined as
// 
// .. math::
//
//     W(\nabla\mathbf{u}) = 
//     \frac{1}{2}\nabla\mathbf{u}\cdot\mathbb{C}\,\nabla\mathbf{u}
//
// For an isotropic material, stress and strain are related through
// 
// .. math::
// 
//     \mathbb{C}_{ijkl} = \lambda \delta_{ij}\varepsilon_{kk} + 2\mu\varepsilon_{ij}
// 
// where :math:`\lambda` and :math:`\mu` are the Lame constant and shear modulus, respectively.
// Users can specify these either through (:code:`lame` and :code:`shear`)
// OR (:code:`lambda` and :code:`mu`) OR (:code:`E` and :code:`nu`).
//
// Class methods:
// 
// #. :code:`Isotropic()`: 
//    Basic constructor. Does nothing, and leaves all values initiated as NAN.
// #. :code:`Isotropic(Solid<Set::Sym::Isotropic> base)`
//    Basic constructor. Does nothing gut allows for inheritance.
// #. :code:`Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda)`
//    BAD old-fashioned constructor. Do not use!
// #. :code:`~Isotropic()`
//    Simple destructor. Don't need to change it.
// #. :code:`void Define(Set::Scalar a_mu, Set::Scalar a_lambda)`
//    BAD old-fashioned way of doing things. Use :code:`Parse` instead.
// #. :code:`Set::Scalar W(const Set::Matrix & gradu) const override`
//    Returns elastic free energy density
// #. :code:`Set::Matrix DW(const Set::Matrix & gradu) const override`
//    Returns first derivative of free energy, the stress tensor
// #. :code:`Set::Matrix4<...> DDW(const Set::Matrix & ) const override`
//    Returns second derivative of free energy, the modulus tensor
// #. :code:`virtual void Print(std::ostream &out) const override`
//    Prints the modulus tensor object to output stream (usually the terminal)
// #. :code:`static Isotropic Random()`
//    Static method that generates a random yet acceptable model.
// #. :code:`static Isotropic Zero()`
//    Static method that generates a "zero" element (so that there is no effect under addition)
// #. :code:`static void Parse(Isotropic & value, IO::ParmParse & pp)`
//    Parser where all the IO occurs
//
//

#ifndef MODEL_SOLID_LINEAR_ISOTROPIC_H_
#define MODEL_SOLID_LINEAR_ISOTROPIC_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Isotropic : public Solid<Set::Sym::Isotropic>
{
public:

    Isotropic() {};
    Isotropic(Solid<Set::Sym::Isotropic> base) : Solid<Set::Sym::Isotropic>(base) {};
    Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda) 
    {
        Define(a_mu,a_lambda);
    };
    virtual ~Isotropic() {};

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
    {
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>(a_lambda,a_mu);
    }

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return ddw*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return ddw;
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << ddw;
    }

    Set::Scalar W2(const Set::Matrix & stress)
    {
        return ( 0.5 * stress.transpose() * (ddw_inv * stress) ).trace();
    }
    Set::Matrix DW2(const Set::Matrix & stress) 
    {
        return ddw_inv * stress;
    }
    Set::Scalar Wvol(const Set::Scalar eps_tr)
    {
        return 0.5*ddw.Lambda()*eps_tr*eps_tr;
    }
    Set::Scalar Wdev(const Set::Matrix eps)
    {
        return ddw.Mu()*std::max((eps*eps).trace(),0.0);
    }

public:
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> ddw;
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> ddw_inv;
    static const KinematicVariable kinvar = KinematicVariable::gradu;

public:
    static Isotropic Random()
    {
        Isotropic ret;
        ret.Define(Util::Random(),Util::Random());
        return ret;
    }
    static Isotropic Zero()
    {
        Isotropic ret;
        ret.Define(0.,0.);
        return ret;
    }
    static void Parse(Isotropic & value, IO::ParmParse & pp)
    {
        Set::Scalar mu = NAN, lambda = NAN;
        bool planestress = false;
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
        
        // Whether or not to use the
        // `plane stress <https://en.wikipedia.org/wiki/Plane_stress>`_ 
        // approximation.
        pp.query_default("planestress",planestress, false); 

        if (AMREX_SPACEDIM==2 && planestress)
            value.Define(mu,lambda*(1.0 - lambda/(2.*mu + lambda)));
        else
            value.Define(mu,lambda);
    }

    #define OP_CLASS Isotropic
    #define OP_VARS X(ddw)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"
}
}
}

template<>
ALAMO_SINGLE_DEFINITION
int Set::Field<Model::Solid::Linear::Isotropic>::NComp() const
{
    return 2;
}
template<>
ALAMO_SINGLE_DEFINITION
std::string Set::Field<Model::Solid::Linear::Isotropic>::Name(int i) const
{
    if (i == 0) return name + "_lambda";
    if (i == 1) return name + "_mu";
    return name;
}
template<>
ALAMO_SINGLE_DEFINITION
void Set::Field<Model::Solid::Linear::Isotropic>::Copy(int a_lev, amrex::MultiFab& a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Linear::Isotropic> const& src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const& dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM * AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i, j, k, a_dstcomp + 0) = src(i, j, k).ddw.Lambda();
                    dst(i, j, k, a_dstcomp + 1) = src(i, j, k).ddw.Mu();
                });
            }
}
}
}

#endif
