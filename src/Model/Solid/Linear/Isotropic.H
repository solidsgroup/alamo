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
    Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda) 
    {
        Define(a_mu,a_lambda);
    };
    virtual ~Isotropic() {};

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
    {
        m_lambda = a_lambda;
        m_mu = a_mu;
        // m_lambda0 = a_lambda;
        // m_mu0 = a_mu;
    }

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (DDW(gradu)*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return DDW(gradu)*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>(m_lambda,m_mu);
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << "lamda = " << m_lambda <<", mu = " << m_mu << std::endl;
        // out << "lamda0 = " << m_lambda0 <<", mu0 = " << m_mu0 << std::endl;
    }

    void DegradeModulus(const Set::Scalar eta)
    {
        if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta ", eta);
        // if(std::isnan(m_lambda0) || std::isnan(m_mu0)) 
        //     Util::Abort(INFO, "Nans detected"," mu0=", m_mu0, " lambda0=",m_lambda0); 

        // m_mu = m_mu0*(1.-eta);
        // m_lambda = m_lambda0*(1.-eta);
    }

    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
    {
        if(eta.size() == 0) return;
        else if(eta.size() == 1) DegradeModulus(eta[0]);
        else if(eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            // m_mu = m_mu0*(1.-eta[0]);
            // m_lambda = m_lambda0*(1.-eta[1]);
            // Isotropic::Define( m_mu, m_lambda );
        }
        else
            return;
    }

public:
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> ddw;
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    Set::Scalar m_lambda = NAN;
    Set::Scalar m_mu = NAN;
    // Set::Scalar m_lambda0 = NAN;
    // Set::Scalar m_mu0 = NAN;

public:
    AMREX_FORCE_INLINE
    static Isotropic Combine(const std::vector<Isotropic> &models, const std::vector<Set::Scalar> &eta)
    {
        Isotropic ret = Isotropic::Zero();
        Set::Scalar etasum = 0.;
        for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        for (unsigned int n = 0 ; n < models.size(); n++)
        {
            ret += models[n] * (eta[n] / etasum);
        }
        return ret;
    }
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
        Set::Scalar mu, lambda;
        bool planestress = false;
        
        // Whether or not to use the
        // `plane stress <https://en.wikipedia.org/wiki/Plane_stress>`_ 
        // approximation.
        pp.query("planestress",planestress); 
        if (pp.contains("lame") && pp.contains("shear"))
        {
            pp.query("lame",lambda); // Lame parameter
            pp.query("shear",mu); // Shear modulus (redundant with "mu")
        }
        if (pp.contains("lambda") && pp.contains("mu"))
        {
            pp.query("lambda",lambda); // Lame parameter
            pp.query("mu",mu); // Shear modulus (redundant with "shear")
        }
        else if (pp.contains("E") && pp.contains("nu"))
        {
            Set::Scalar E, nu;
            pp.query("E",E); // Elastic modulus
            pp.query("nu",nu); // Poisson's ratio
            lambda = E * nu / (1.0 + nu) / (1.0 - 2.0*nu);
            mu = E / 2.0 / (1.0 + nu);
        }
        else Util::Abort(INFO,"Invalid elastic constants specified");

        if (AMREX_SPACEDIM==2 && planestress)
            value.Define(mu,lambda*(1.0 - lambda/(2.*mu + lambda)));
        else
            value.Define(mu,lambda);
        Util::Message(INFO,value);
    }

    #define OP_CLASS Isotropic
    #define OP_VARS X(m_lambda) X(m_mu) //X(m_lambda0) X(m_mu0)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"


}
}
}

template<>
inline int Set::Field<Model::Solid::Linear::Isotropic>::NComp() const
{
    return 2;    
}
template<>
inline std::string Set::Field<Model::Solid::Linear::Isotropic>::Name(int i) const
{
    if (i==0) return name + "_lambda";
    if (i==1) return name + "_mu";
    return name;
}
template<>
inline void Set::Field<Model::Solid::Linear::Isotropic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Linear::Isotropic> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM*AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + 0) = src(i,j,k).m_lambda;
                    dst(i,j,k,a_dstcomp + 1) = src(i,j,k).m_mu;
                });
            }
        }
    }    
}

#endif
