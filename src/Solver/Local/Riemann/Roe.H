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

/// Roe Riemann Solver based on Gas Dynamics - Culbert B. Laney
class Roe
{
public:
    //std::array<double, 4> Roe(double left_state[5], double right_state[5], Set::Scalar gamma) {
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
        Set::Scalar u_lo = lo.momentum(0) / lo.rho;
        Set::Scalar ke_lo = 0.5 * lo.rho * u_lo * u_lo;
        Set::Scalar ue_lo = lo.energy - ke_lo;
        Set::Scalar p_lo = (gamma - 1) * ue_lo;
        Set::Scalar h_lo = (lo.energy + p_lo) / lo.rho;

        Set::Scalar u_hi = hi.momentum(0) / hi.rho;
        Set::Scalar ke_hi = 0.5 * hi.rho * u_hi * u_hi;
        Set::Scalar ue_hi = hi.energy - ke_hi;
        Set::Scalar p_hi = (gamma - 1) * ue_hi;
        Set::Scalar h_hi = (hi.energy + p_hi) / hi.rho;

        // roe averages
        Set::Scalar rho_RL = std::sqrt(lo.rho) * std::sqrt(hi.rho);
        Set::Scalar u_RL   = (std::sqrt(hi.rho) * u_hi + std::sqrt(lo.rho) * u_lo) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));
	Set::Scalar h_RL   = (std::sqrt(hi.rho) * h_hi + std::sqrt(lo.rho) * h_lo) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));
        Set::Scalar a_RL2  = (gamma - 1) * (h_RL - 0.5 * u_RL * u_RL);
        Set::Scalar a_RL   = std::sqrt(a_RL2);
        Set::Scalar eta_RL = (std::sqrt(hi.rho) * hi.eta + std::sqrt(lo.rho) * lo.eta) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));

        // flux
        Flux fl;
        fl.mass =        0.5 * (lo.rho * u_lo * lo.eta + hi.rho * u_hi * hi.eta);
        fl.momentum(0) = 0.5 * (lo.rho * u_lo * u_lo * lo.eta + hi.rho * u_hi * u_hi * hi.eta) + 0.5 * (p_lo + p_hi) * eta_RL;
        fl.momentum(1) = 0.0;
        fl.energy =      0.5 * (u_lo * ke_lo * lo.eta + u_hi * ke_hi * hi.eta) + 0.5 * (u_lo * (ue_lo + p_lo) + u_hi * (ue_hi + p_hi)) * eta_RL;

	//compute eigenvalues
	std::array <double, 3> ev;
        ev[0] = std::abs(u_RL);
        ev[1] = std::abs(u_RL + a_RL);
        ev[2] = std::abs(u_RL - a_RL);
	
	//compute wave strengths

	Set::Scalar drho, dp, dvel;
	if (p_hi < p_lo) {
	dvel = u_hi * hi.eta - u_lo * lo.eta;
	dp   = p_hi - p_lo;
	drho = hi.rho * hi.eta - lo.rho * lo.eta;
	}
	else {
	dvel = -(u_hi * hi.eta - u_lo * lo.eta);
	dp   = -(p_hi - p_lo);
	drho = -(hi.rho * hi.eta - lo.rho * lo.eta);
	};

	std::array <double, 3> dv;
        dv[0] = drho - dp / a_RL2;
        dv[1] = dvel + dp / (rho_RL * a_RL);
        dv[2] = dvel - dp / (rho_RL * a_RL);
	
	//construct right characteristic eigenvectors
        std::array <double, 3> rev_0;
        rev_0[0] = 1;
	rev_0[1] = u_RL;
	rev_0[2] = 0.5 * u_RL * u_RL;

	std::array <double, 3> rev_1;
	rev_1[0] = (rho_RL / 2 / a_RL);
	rev_1[1] = (rho_RL / 2 / a_RL) * (u_RL + a_RL);
	rev_1[2] = (rho_RL / 2 / a_RL) * (h_RL + a_RL * u_RL);

	std::array <double, 3> rev_2;
	rev_2[0] = -(rho_RL / 2 / a_RL);
	rev_2[1] = -(rho_RL / 2 / a_RL) * (u_RL - a_RL);
	rev_2[2] = -(rho_RL / 2 / a_RL) * (h_RL - a_RL * u_RL);
	
	//upwinding term 
	fl.mass        += -0.5 * (rev_0[0] * ev[0] * dv[0] + rev_1[0] * ev[1] * dv[1] + rev_2[0] * ev[2] * dv[2]);
	fl.momentum(0) += -0.5 * (rev_0[1] * ev[0] * dv[0] + rev_1[1] * ev[1] * dv[1] + rev_2[1] * ev[2] * dv[2]);
	fl.energy      += -0.5 * (rev_0[2] * ev[0] * dv[0] + rev_1[2] * ev[1] * dv[1] + rev_2[2] * ev[2] * dv[2]);

        return fl;
    }

    static int Test()
    {
        int failed = 0;
        for (int i = 0; i < 10; i++)
        {

	    State hi, lo, fixed, pos_vel1, pos_vel2, neg_vel1, neg_vel2;
            hi.rho = Util::Random();
            hi.momentum = Set::Vector::Random();
            hi.energy = Util::Random();
            hi.eta = Util::Random();

            lo.rho = Util::Random();
            lo.momentum = Set::Vector::Random();
            lo.energy = Util::Random();
            lo.eta = Util::Random();

            lo.rho = Util::Random();
            lo.momentum = Set::Vector::Zero();
            lo.energy = Util::Random();
            lo.eta = Util::Random();

	    pos_vel1.rho = Util::Random();
	    pos_vel1.momentum = Set::Vector::Random();
	    pos_vel1.energy = Util::Random();
	    pos_vel1.eta = Util::Random();

	    neg_vel1.rho = pos_vel1.rho;
	    neg_vel1.momentum = -pos_vel1.momentum;
	    neg_vel1.energy = pos_vel1.energy;
	    neg_vel1.eta = pos_vel1.eta;

	    pos_vel2.rho = Util::Random();
	    pos_vel2.momentum = Set::Vector::Random();
	    pos_vel2.energy = Util::Random();
	    pos_vel2.eta = Util::Random();

	    neg_vel2.rho = pos_vel2.rho;
	    neg_vel2.momentum = -pos_vel2.momentum;
	    neg_vel2.energy = pos_vel2.energy;
	    neg_vel2.eta = pos_vel2.eta;

            // 
            // Zero flux test: if the state is the same on both sides,
            //                 the flux should be zero.
            // 
            Flux zero = Solve(fixed,fixed,1.4);
            if (fabs(zero.mass) > 1E-10)  {Util::Warning(INFO,"Nonzero mass flux ", zero.mass); failed++;}
            if (fabs(zero.momentum(0)) > 1E-10)  {Util::Warning(INFO,"Nonzero x momentum flux ", zero.momentum(0)); failed++;}
            if (fabs(zero.momentum(1)) > 1E-10)  {Util::Warning(INFO,"Nonzero y momentum flux ", zero.momentum(1)); failed++;}
            if (fabs(zero.energy) > 1E-10)  {Util::Warning(INFO,"Nonzero energy flux ", zero.energy); failed++;}

            // 
            // Antisymmetry test: reversing hi and lo should produce the same
            //                    output 
            // 

            Flux plus = Solve(hi, lo, 1.4);
            Flux minus = Solve(lo, hi, 1.4);

            if (fabs(plus.mass - minus.mass) > 1E-10)  
            {Util::Warning(INFO,"mass flux error: plus=", plus.mass, " minus=",minus.mass); failed++;}
            if (fabs(plus.momentum(0) - minus.momentum(0)) > 1E-10)  
            {Util::Warning(INFO,"momentum x flux error: plus=", plus.momentum(0), " minus=",minus.momentum(0)); failed++;}
            if (fabs(plus.momentum(1) - minus.momentum(1)) > 1E-10)  
            {Util::Warning(INFO,"momentum y flux error: plus=", plus.momentum(1), " minus=",minus.momentum(1)); failed++;}
            if (fabs(plus.energy - minus.energy) > 1E-10)  
            {Util::Warning(INFO,"energy flux error: plus=", plus.energy, " minus=",minus.energy); failed++;}


	    // 
            // Wave Propogation Test: changing the sign of velocity should produce negative mass
            //                   and energy fluxes and the same momentum flux
            // 

            Flux pos = Solve(pos_vel1, pos_vel2, 1.4);
            Flux neg = Solve(neg_vel1, neg_vel2, 1.4);

            if (fabs(pos.mass + neg.mass) > 1E-10)  
            {Util::Warning(INFO,"mass flux error: positive vel=", pos.mass, " negative vel=", neg.mass); failed++;}
            if (fabs(pos.momentum(0) - neg.momentum(0)) > 1E-10)  
            {Util::Warning(INFO,"momentum x flux error: positive vel=", pos.momentum(0), " negative vel=", neg.momentum(0)); failed++;}
            if (fabs(pos.momentum(1) - neg.momentum(1)) > 1E-10)  
            {Util::Warning(INFO,"momentum y flux error: positive vel=", pos.momentum(1), " negative vel=", neg.momentum(1)); failed++;}
            if (fabs(pos.energy + neg.energy) > 1E-10)  
            {Util::Warning(INFO,"energy flux error: positive vel=", pos.energy, " negative vel=", neg.energy); failed++;}


        }

        return failed;
    }
};
}
}
}
