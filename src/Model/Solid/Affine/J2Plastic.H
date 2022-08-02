#ifndef MODEL_SOLID_PLASTIC_J2_H_
#define MODEL_SOLID_PLASTIC_J2_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class J2Plastic : public Isotropic
{
public:
    J2Plastic() {};
    J2Plastic(Solid<Set::Sym::Isotropic> base) : Isotropic(base) {};

    J2Plastic(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta=1.0)
    {
        Define(a_mu, a_lambda, a_yield, a_hardening, a_theta);
    }

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta)
    {
        m_mu = a_mu;
        theta = a_theta;
        yield_strength = a_yield;
        hardening_modulus = a_hardening;
        
        curr = PlasticState::Zero(); 
        prev = PlasticState::Zero();

        Isotropic::Define(a_mu, a_lambda);
    }

    Set::Scalar YieldSurface()
    {
        return yield_strength + curr.alpha*hardening_modulus;
    }

    Set::Scalar PlasticEnergy()
    {
        return (yield_strength*curr.alpha + 0.5*hardening_modulus*curr.alpha*curr.alpha);
    }
    
    void EvolvePlasticStrain(Set::Matrix sigma, Set::Matrix strain, Set::Scalar /*dt*/)
    {
        Set::Scalar SQ2O3 = sqrt(1.0 - 1.0/((double)AMREX_SPACEDIM));
        Set::Matrix sigdev = sigma - (1.0/((double)AMREX_SPACEDIM))*sigma.trace()*Set::Matrix::Identity();
        Set::Matrix epsdev = strain - strain.trace()*Set::Matrix::Identity();

        Set::Matrix zeta_trial = sigdev  - prev.beta + 2.0*m_mu*epsdev;

        Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(yield_strength + theta*hardening_modulus*prev.alpha);
        if( f_trial <= 0.0)
        {
            curr.epsp = prev.epsp; curr.alpha = prev.alpha; curr.beta = prev.beta; 
            return;
        }
        Set::Matrix n_new = zeta_trial/zeta_trial.norm();
        Set::Scalar dGamma = f_trial/(2.0*(m_mu)*(1.0 + (hardening_modulus/(3.0*m_mu))));
        Set::Scalar dH = SQ2O3*(1.0-theta)*hardening_modulus*dGamma;

        curr.alpha = prev.alpha + SQ2O3*dGamma;
        curr.beta = prev.beta + SQ2O3*dH*n_new;
        curr.epsp = prev.epsp + dGamma*n_new;
        Isotropic::SetF0(curr.epsp);
    }

    void SetPlasticStrains()
    {
        prev = Affine::PlasticState::Zero();
        curr = Affine::PlasticState::Zero();
    }

    void SetPlasticStrains(PlasticState &a_state)
    {
        prev = a_state;
        curr = a_state;
        Isotropic::SetF0(curr.epsp);
    }

    PlasticState GetPlasticState()
    {
        return curr;
    }

public:
    PlasticState curr, prev;
    Set::Scalar theta;                  // isotropic and kinematic hardening parameter
    Set::Scalar yield_strength;          // yield strength
    Set::Scalar hardening_modulus;      // hardening modulus
    Set::Scalar m_mu;

public:
    static void Parse(J2Plastic & value, IO::ParmParse & pp)
    {
        Set::Scalar mu, lambda;
        if (pp.contains("lambda") && pp.contains("mu"))
        {
            pp.query("lambda",lambda);
            pp.query("mu",mu);
        }
        else if (pp.contains("E") && pp.contains("nu"))
        {
            Set::Scalar E, nu;
            pp.query("E",E);
            pp.query("nu",nu);
            lambda = E * nu / (1.0 + nu) / (1.0 - 2.0*nu);
            mu = E / 2.0 / (1.0 + nu);
        }
        Set::Scalar yield = 1.0, hardening = 1.0, theta_tmp = 1.0;

        pp.query("yield", yield);
        pp.query("hardening", hardening);
        pp.query("theta", theta_tmp);

        value.Define(mu, lambda, yield, hardening, theta_tmp);
    }
    AMREX_FORCE_INLINE
    void operator += (const J2Plastic &rhs)
    {
        m_mu += rhs.m_mu;
        curr += rhs.curr;
        prev += rhs.prev;
        theta += rhs.theta;
        yield_strength += rhs.yield_strength;
        hardening_modulus += rhs.hardening_modulus;
        
        ddw += rhs.ddw;
        F0 += rhs.F0;
    }    
    AMREX_FORCE_INLINE
        J2Plastic operator * (const Set::Scalar alpha) const
    {
        J2Plastic ret;
        ret.m_mu = alpha * m_mu;
        ret.curr = alpha*curr;
        ret.prev = alpha*prev;
        ret.theta = alpha*theta;
        ret.yield_strength = alpha*yield_strength;
        ret.hardening_modulus = alpha*hardening_modulus;

        ret.ddw = ddw*alpha;
        ret.F0 = F0*alpha;
    
        return ret;
    }    
    friend J2Plastic operator * (const Set::Scalar alpha, const J2Plastic b);
    friend J2Plastic operator + (const J2Plastic a, const J2Plastic b);
    friend J2Plastic operator - (const J2Plastic a, const J2Plastic b);

};
AMREX_FORCE_INLINE
J2Plastic operator * (const Set::Scalar alpha, const J2Plastic b)
{
    J2Plastic ret;
    ret.m_mu = alpha * b.m_mu;
    ret.curr = alpha * b.curr;
    ret.prev = alpha * b.prev;
    ret.theta = alpha * b.theta;
    ret.yield_strength = alpha * b.yield_strength;
    ret.hardening_modulus = alpha * b.hardening_modulus;

    ret.ddw = b.ddw*alpha;
    ret.F0 = b.F0*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
J2Plastic operator + (const J2Plastic a, const J2Plastic b)
{
    J2Plastic ret;
    ret.m_mu = a.m_mu + b.m_mu;
    ret.curr = a.curr + b.curr;
    ret.prev = a.prev + b.prev;
    ret.theta = a.theta + b.theta;
    ret.yield_strength = a.yield_strength + b.yield_strength;
    ret.hardening_modulus = a.hardening_modulus + b.hardening_modulus;
    
    ret.ddw = a.ddw + b.ddw;
    ret.F0 = a.F0 + b.F0;
    return ret;
}    
AMREX_FORCE_INLINE
J2Plastic operator - (const J2Plastic a, const J2Plastic b)
{
    J2Plastic ret;
    ret.m_mu = a.m_mu - b.m_mu;
    ret.curr = a.curr - b.curr;
    ret.prev = a.prev - b.prev;
    ret.theta = a.theta - b.theta;
    ret.yield_strength = a.yield_strength - b.yield_strength;
    ret.hardening_modulus = a.hardening_modulus - b.hardening_modulus;
    
    ret.ddw = a.ddw - b.ddw;
    ret.F0 = a.F0 - b.F0;
    return ret;
} 


}
}
}
#endif
