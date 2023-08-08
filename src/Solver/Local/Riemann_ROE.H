#include <eigen3/Eigen/Core>
#include "Set/Set.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{
/// Roe Riemann Solver based on Gas Dynamics - Culbert B. Laney
  std::array<double, 4> Riemann_ROE(double left_state[5], double right_state[5], Set::Scalar gamma) {
    Set::Scalar rho_l = left_state[0];
    Set::Scalar u_l = left_state[1]/left_state[0];
    Set::Scalar e_l = left_state[3];
    Set::Scalar ke_l = 0.5 * rho_l * u_l * u_l;
    Set::Scalar ue_l = e_l - ke_l;
    Set::Scalar p_l = (gamma - 1) * ue_l;
    Set::Scalar eta_l = left_state[4];
    Set::Scalar h_l = (e_l + p_l) / rho_l;

    Set::Scalar rho_r = right_state[0];
    Set::Scalar u_r = right_state[1]/right_state[0];
    Set::Scalar e_r = right_state[3];
    Set::Scalar ke_r = 0.5 * rho_r * u_r * u_r;
    Set::Scalar ue_r = e_r - ke_r;
    Set::Scalar p_r = (gamma - 1) * ue_r;
    Set::Scalar eta_r = right_state[4];
    Set::Scalar h_r = (e_r + p_r) / rho_r;

    //Util::Message()

    // roe averages
    Set::Scalar srl = std::sqrt(rho_l);
    Set::Scalar srr = std::sqrt(rho_r);
    Set::Scalar rho_RL = srl * srr;
    Set::Scalar u_RL = (srr * u_r + srl * u_l) / (srl + srr);
    Set::Scalar h_RL = (srr * h_r + srl * h_l) / (srl + srr);
    Set::Scalar a_RL2 = (gamma - 1) * (h_RL - 0.5 * u_RL * u_RL);
    Set::Scalar a_RL = std::sqrt(a_RL2);
    Set::Scalar eta_RL = (srr * eta_r + srl * eta_l) / (srl + srr);

    // eigenvalues
    std::array<Set::Scalar, 4> ev;
    ev[0] = u_RL;
    ev[2] = u_RL + a_RL;
    ev[3] = 0;
    ev[1] = u_RL - a_RL;
    
    std::array<Set::Scalar, 4> dv;
    Set::Scalar drho = rho_r - rho_l;
    Set::Scalar dp = p_r - p_l;
    Set::Scalar dvel = u_r - u_l;
    dv[0] = drho - dp / a_RL2;
    dv[2] = dvel + dp / (rho_RL * a_RL);
    dv[3] = 0;
    dv[1] = dvel - dp / (rho_RL * a_RL);

    // right_state characteristic eigenvectors
    std::array<Set::Scalar, 4> rev_0;
    rev_0[0] = 1;
    rev_0[2] = u_RL;
    rev_0[3] = 0;
    rev_0[1] = 0.5 * u_RL * u_RL;

    std::array<Set::Scalar, 4> rev_1;
    rev_1[0] = (rho_RL / 2. / a_RL);
    rev_1[2] = (rho_RL / 2. / a_RL) * (u_RL + a_RL);
    rev_1[3] = 0;
    rev_1[1] = (rho_RL / 2. / a_RL) * (h_RL + a_RL * u_RL);

    std::array<Set::Scalar, 4> rev_2;
    rev_2[0] = 0;
    rev_2[2] = 0;
    rev_2[3] = 0;
    rev_2[1] = 0;

    std::array<Set::Scalar, 4> rev_3;
    rev_3[0] = (-rho_RL / 2 / a_RL);
    rev_3[2] = (-rho_RL / 2 / a_RL) * (u_RL - a_RL);
    rev_3[3] = 0;
    rev_3[1] = (-rho_RL / 2 / a_RL) * (h_RL - a_RL * u_RL);

    // flux
    std::array<Set::Scalar, 4> fl;
    fl[0] = 0.5 * (rho_l * u_l * eta_l + rho_r * u_r * eta_r);
    fl[2] = 0.5 * (rho_l * u_l * u_l * eta_l + rho_r * u_r * u_r * eta_r) + 0.5 * (p_l + p_r) * eta_RL;
    fl[3] = 0.0;
    fl[1] = 0.5 * (u_l * ke_l * eta_l + u_r * ke_r * eta_r) + 0.5 * (u_l * (ue_l + p_l) + u_r * (ue_r + p_r)) * eta_RL;

    fl[0] += -0.5 * rev_0[0] * std::abs(ev[0]) * dv [0];
    fl[0] += -0.5 * rev_1[0] * std::abs(ev[0]) * dv [0];
    fl[0] += -0.5 * rev_2[0] * std::abs(ev[0]) * dv [0];
    fl[0] += -0.5 * rev_3[0] * std::abs(ev[0]) * dv [0];
    
    fl[2] += -0.5 * rev_0[2] * std::abs(ev[2]) * dv [2];
    fl[2] += -0.5 * rev_1[2] * std::abs(ev[2]) * dv [2];
    fl[2] += -0.5 * rev_2[2] * std::abs(ev[2]) * dv [2];
    fl[2] += -0.5 * rev_3[2] * std::abs(ev[2]) * dv [2];
    
    fl[1] += -0.5 * rev_0[1] * std::abs(ev[1]) * dv [1];
    fl[1] += -0.5 * rev_1[1] * std::abs(ev[1]) * dv [1];
    fl[1] += -0.5 * rev_2[1] * std::abs(ev[1]) * dv [1];
    fl[1] += -0.5 * rev_3[1] * std::abs(ev[1]) * dv [1];
    
    return fl;
}
}
}
