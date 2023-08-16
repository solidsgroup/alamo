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
    static Flux Solve(State left, State right, Set::Scalar gamma) {
        Set::Scalar rho_l = left.rho;
        Set::Scalar u_l = left.momentum(0) / left.rho;
        Set::Scalar e_l = left.energy;
        Set::Scalar ke_l = 0.5 * rho_l * u_l * u_l;
        Set::Scalar ue_l = e_l - ke_l;
        Set::Scalar p_l = (gamma - 1) * ue_l;
        Set::Scalar eta_l = left.eta;
        Set::Scalar h_l = (e_l + p_l) / rho_l;

        Set::Scalar rho_r = right.rho;
        Set::Scalar u_r = right.momentum(0) / right.rho;
        Set::Scalar e_r = right.energy;
        Set::Scalar ke_r = 0.5 * rho_r * u_r * u_r;
        Set::Scalar ue_r = e_r - ke_r;
        Set::Scalar p_r = (gamma - 1) * ue_r;
        Set::Scalar eta_r = right.eta;
        Set::Scalar h_r = (e_r + p_r) / rho_r;

        // roe averages
        Set::Scalar srl = std::sqrt(rho_l);
        Set::Scalar srr = std::sqrt(rho_r);
        Set::Scalar rho_RL = srl * srr;
        Set::Scalar u_RL = (srr * u_r + srl * u_l) / (srl + srr);
        Set::Scalar h_RL = (srr * h_r + srl * h_l) / (srl + srr);
        Set::Scalar a_RL2 = (gamma - 1) * (h_RL - 0.5 * u_RL * u_RL);
        Set::Scalar a_RL = std::sqrt(a_RL2);
        Set::Scalar eta_RL = (srr * eta_r + srl * eta_l) / (srl + srr);

	Set::Scalar sgn;
	if (u_RL > 0) {sgn = 1;}
	if (u_RL < 0) {sgn = -1;}
	else {sgn = 0;};

        // eigenvalues
        std::array<Set::Scalar, 4> ev;
        ev[0] = u_RL;
        ev[2] = u_RL + a_RL;
        ev[3] = 0;
        ev[1] = u_RL - a_RL;

        std::array<Set::Scalar, 4> dv;

        Set::Scalar drho = sgn * std::abs(rho_r * eta_r - rho_l * eta_l);
        Set::Scalar dp   = sgn * std::abs(p_r - p_l);
        Set::Scalar dvel = sgn * std::abs(u_r * eta_r - u_l * eta_l);

        dv[0] = drho - dp * eta_RL / a_RL2;
        dv[2] = dvel + dp * eta_RL / (rho_RL * a_RL);
        dv[3] = 0;
        dv[1] = dvel - dp * eta_RL / (rho_RL * a_RL);

        // right_state characteristic eigenvectors
        std::array<Set::Scalar, 4> rev_0;
        rev_0[0] = 1;
        rev_0[2] = u_RL;
        rev_0[3] = 0;
        rev_0[1] = 0.5 * u_RL * u_RL;

        std::array<Set::Scalar, 4> rev_1;
        rev_1[0] = (rho_RL / 2. / a_RL) * eta_RL;
        rev_1[2] = (rho_RL / 2. / a_RL) * (u_RL + a_RL) * eta_RL;
        rev_1[3] = 0;
        rev_1[1] = (rho_RL / 2. / a_RL) * (h_RL + a_RL * u_RL) * eta_RL;

        std::array<Set::Scalar, 4> rev_2;
        rev_2[0] = 0;
        rev_2[2] = 0;
        rev_2[3] = 0;
        rev_2[1] = 0;

        std::array<Set::Scalar, 4> rev_3;
        rev_3[0] = (-rho_RL / 2 / a_RL) * eta_RL;
        rev_3[2] = (-rho_RL / 2 / a_RL) * (u_RL - a_RL) * eta_RL;
        rev_3[3] = 0;
        rev_3[1] = (-rho_RL / 2 / a_RL) * (h_RL - a_RL * u_RL) * eta_RL;

        // flux
        //std::array<Set::Scalar, 4> fl;
        Flux fl;
        fl.mass = 0.5 * (rho_l * u_l * eta_l + rho_r * u_r * eta_r);
        fl.momentum(0) = 0.5 * (rho_l * u_l * u_l * eta_l + rho_r * u_r * u_r * eta_r) + 0.5 * (p_l + p_r) * eta_RL;
        fl.momentum(1) = 0.0;
        fl.energy = 0.5 * (u_l * ke_l * eta_l + u_r * ke_r * eta_r) + 0.5 * (u_l * (ue_l + p_l) + u_r * (ue_r + p_r)) * eta_RL;

        // upwinding

        // fl.mass += -0.5 * rev_0[0] * std::abs(ev[0]) * dv[0];
        // fl.mass += -0.5 * rev_1[0] * std::abs(ev[1]) * dv[1];
        // fl.mass += -0.5 * rev_2[0] * std::abs(ev[2]) * dv[2];
        // fl.mass += -0.5 * rev_3[0] * std::abs(ev[3]) * dv[3];

        // fl.momentum(0) += -0.5 * rev_0[2] * std::abs(ev[0]) * dv[0];
        // fl.momentum(0) += -0.5 * rev_1[2] * std::abs(ev[1]) * dv[1];
        // fl.momentum(0) += -0.5 * rev_2[2] * std::abs(ev[2]) * dv[2];
        // fl.momentum(0) += -0.5 * rev_3[2] * std::abs(ev[3]) * dv[3];

        // fl.energy += -0.5 * rev_0[1] * std::abs(ev[0]) * dv[0];
        // fl.energy += -0.5 * rev_1[1] * std::abs(ev[1]) * dv[1];
        // fl.energy += -0.5 * rev_2[1] * std::abs(ev[2]) * dv[2];
        // fl.energy += -0.5 * rev_3[1] * std::abs(ev[3]) * dv[3];

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
