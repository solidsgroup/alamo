#include <eigen3/Eigen/Core>
#include "Set/Set.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{

namespace Riemann
{

/// HLLC Riemann Solver
class HLLC
{
public:
    struct State {
        Set::Scalar rho;
        Set::Vector momentum;
        Set::Scalar energy;
        Set::Scalar eta;
        // Construtor for convenience
        State() {rho=0.0; momentum=Set::Vector::Zero(); energy = 0.0; eta = 0.0;} 
        State(Set::Scalar a_rho, Set::Scalar a_mx, Set::Scalar a_my, Set::Scalar a_energy, Set::Scalar a_eta)
            : rho(a_rho), momentum(a_mx, a_my), energy(a_energy), eta(a_eta) {}
    };

    struct Flux {
        Set::Scalar mass;
        Set::Scalar energy;
        Set::Vector momentum;
    };


    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma) {
      
	//lo primitive variables
        Set::Scalar uNorm_lo = lo.momentum(0) / lo.rho;
	Set::Scalar uTang_lo = lo.momentum(1) / lo.rho;
        Set::Scalar ke_lo = 0.5 * lo.rho * (uNorm_lo * uNorm_lo + uTang_lo * uTang_lo);
        Set::Scalar ue_lo = lo.energy - ke_lo;
        Set::Scalar p_lo = (gamma - 1) * ue_lo;

	//hi primitive variables
        Set::Scalar uNorm_hi = hi.momentum(0) / hi.rho;
	Set::Scalar uTang_hi = hi.momentum(1) / hi.rho;
        Set::Scalar ke_hi = 0.5 * lo.rho * (uNorm_hi * uNorm_hi + uTang_hi * uTang_hi);
        Set::Scalar ue_hi = hi.energy - ke_hi;
        Set::Scalar p_hi = (gamma - 1) * ue_hi;

	//Find the largest eigenvalues in the normal direction to the interface
        Set::Scalar a_lo = std::sqrt(gamma * p_lo/lo.rho);
	Set::Scalar a_hi = std::sqrt(gamma * p_hi/hi.rho);

	//Compute HLL wave speed
        Set::Scalar S_lo = std::min(uNorm_lo, uNorm_hi) - max(a_lo, a_hi);
	Set::Scalar S_hi = std::max(uNorm_lo, uNorm_hi) + max(a_lo, a_hi);

	//Compute lagrangian sound speed
        Set::Scalar aLag_lo = lo.rho * (uNorm_lo - S_lo);
        Set::Scalar aLag_hi = hi.rho * (S_hi - uNorm_hi);

	//Compute acoustic star state
        Set::Scalar ustar = (aLag_hi * uNorm_hi + aLag_lo * uNorm_lo + (p_lo - p_hi))/(aLag_hi + aLag_lo);
        Set::Scalar pstar = (aLag_hi * p_lo + aLag_lo * p_hi + aLag_lo * aLag_hi * (uNorm_lo - uNorm_hi))/(aLag_hi + aLag_lo);

	//Left star region variables
        Set::Scalar rstar_lo   = lo.rho * (S_lo - uNorm_lo)/(S_lo - ustar);
        Set::Scalar estar_lo   = ((S_lo - uNorm_lo) * lo.energy - p_lo * uNorm_lo + pstar * ustar)/(S_lo - ustar);
	Set::Scalar etastar_lo = lo.eta * (S_lo - uNorm_lo)/(S_lo - ustar);

	//Right star region variables
        Set::Scalar rstar_hi   = hi.rho * (S_hi - uNorm_hi)/(S_hi - ustar);
        Set::Scalar estar_hi   = ((S_hi - uNorm_hi) * hi.energy - p_hi * uNorm_hi + pstar * ustar)/(S_hi - ustar);
	Set::Scalar etastar_hi = hi.eta * (S_hi - uNorm_hi)/(S_hi - ustar);

	//Sample the solution at x/t=0
	Set::Scalar rho_face, u_face, p_face, e_face, eta_face;
        if (S_lo > 0.0) {
	  rho_face = lo.rho;
          u_face   = uNorm_lo;
          p_face   = p_lo;
          e_face   = lo.energy;
	  eta_face = lo.eta;
	}
        else if (ustar > 0.0) {
	  rho_face = rstar_lo;
          u_face   = ustar;
          p_face   = pstar;
          e_face   = estar_lo;
	  eta_face = etastar_lo;
	}
        else if (S_hi > 0.0) {
	  rho_face = rstar_hi;
          u_face = ustar;
          p_face = pstar;
          e_face = estar_hi;
	  eta_face = etastar_hi;
        }
        else {
	  rho_face = hi.rho;
          u_face = uNorm_hi;
          p_face = p_hi;
          e_face = hi.energy;
	  eta_face = hi.eta;
	};

	//Compute the Godunov flux
        Flux fl;
        fl.mass = rho_face * u_face * eta_face;
        fl.momentum(0) = (rho_face * u_face * u_face) * eta_face + p_face;
        if (fl.mass > 0.0) {
	  fl.momentum(1) = (fl.mass * uTang_lo) * u_face * eta_face;
	}
	else {
	  fl.momentum(1) = (fl.mass * uTang_hi) * u_face * eta_face;
	};
        fl.energy = (e_face + p_face) * u_face * eta_face;

        return fl;
	      
    }
};
}
}
}
