#ifndef MODEL_SOLID_PLASTIC_J2_H_
#define MODEL_SOLID_PLASTIC_J2_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
template <class T>
class J2Plastic : public T
{
public:
    J2Plastic() {};
    J2Plastic(Solid<Set::Sym::Isotropic> base) : T(base) {};

    J2Plastic(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta=1.0)
    {
        Define(a_mu, a_lambda, a_yield, a_hardening, a_theta);
    }

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta)
    {
        theta = a_theta;
        yield_strength = a_yield;
        yield_strength0 = a_yield;
        hardening_modulus = a_hardening;
        hardening_modulus0 = a_hardening;

        curr.epsp = Set::Matrix::Zero(); 
        prev.epsp = Set::Matrix::Zero();
        curr.beta = Set::Matrix::Zero();
        prev.beta = Set::Matrix::Zero();
        curr.alpha = 0.0; 
        prev.alpha = 0.0;

        T::Define(a_mu, a_lambda);
    }

    void DegradeYieldSurface(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
        if(std::isinf(eta)) Util::Abort(INFO, "Infs in eta");
		yield_strength = yield_strength0*(1. - eta);
        hardening_modulus = hardening_modulus0*(1.-eta);
	}

    void DegradeYieldSurface(const amrex::Vector<Set::Scalar> &eta)
    {
        if (eta.size() == 0) return;
        else if (eta.size() == 1)
        {
            if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
			yield_strength = yield_strength0*(1. - eta[0]);
            hardening_modulus = hardening_modulus0*(1.-eta[0]);
        }
        else if (eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            yield_strength = yield_strength0*(1. - eta[0]);
            hardening_modulus = hardening_modulus0*(1.-eta[1]);
        }
        else return;
    }

    Set::Scalar YieldSurface()
    {
        return yield_strength + curr.alpha*hardening_modulus;
    }

    Set::Scalar OriginalYieldSurface()
    {
        return yield_strength0 + curr.alpha*hardening_modulus0;
    }

    void EvolvePlasticStrain(Set::Matrix sigma, Set::Matrix /*strain*/, Set::Scalar /*dt*/)
    {
        Set::Scalar SQ2O3 = sqrt(1.0 - 1.0/((double)AMREX_SPACEDIM));
        Set::Matrix sigdev = sigma - (1.0/((double)AMREX_SPACEDIM))*sigma.trace()*Set::Matrix::Identity();
        //Set::Matrix epsdev = strain - strain.trace()*Set::Matrix::Identity();

        Set::Matrix zeta_trial = sigdev - prev.beta;

        Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(yield_strength + theta*hardening_modulus*prev.alpha);
        if( f_trial < 0.0)
        {
            curr.epsp = curr.epsp; curr.alpha = prev.alpha; curr.beta = prev.beta; 
            return;
        }
        Set::Matrix n_new = zeta_trial/zeta_trial.norm();
        Set::Scalar dGamma = f_trial/(2.0*(T::mu)*(1.0 + (hardening_modulus/(3.0*T::mu))));
        Set::Scalar dH = SQ2O3*(1.0-theta)*hardening_modulus*dGamma;

        curr.alpha = prev.alpha + SQ2O3*dGamma;
        curr.beta = prev.beta + SQ2O3*dH*n_new;
        curr.epsp = prev.epsp + dGamma*n_new;
        T::UpdateF0(curr.epsp);
    }

    void SetPlasticStrains(PlasticState &a_state)
    {
        prev.epsp = a_state.epsp; curr.epsp = a_state.epsp;
        prev.beta = a_state.beta; curr.beta = a_state.beta;
        prev.alpha = a_state.alpha; curr.alpha = a_state.alpha;
    }

    PlasticState GetPlasticState()
    {
        return curr;
    }

public:
    PlasticState curr, prev;
    Set::Scalar theta;                  // isotropic and kinematic hardening parameter
    Set::Scalar yield_strength, yield_strength0;          // yield strength
    Set::Scalar hardening_modulus, hardening_modulus0;  // hardening modulus

public:
    static void Parse(J2Plastic<T> & value, IO::ParmParse & pp)
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
};

}
}
}
#endif