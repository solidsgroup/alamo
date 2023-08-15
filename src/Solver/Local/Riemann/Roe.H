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
std::array<double, 4> Roe(double left_state[5], double right_state[5], Set::Scalar eta, Set::Scalar gamma) {
    Set::Scalar rho_l, vel_l, p_l, ke_l, ue_l, h_l;
    Set::Scalar rho_r, vel_r, p_r, ke_r, ue_r, h_r;
    Set::Scalar eta_r, eta_l;

    rho_l = left_state[0];
    vel_l = left_state[1];
    p_l = left_state[3];
    eta_l = left_state[4];

    ke_l = 0.5 * rho_l * vel_l * vel_l;
    ue_l = p_l / (gamma - 1.0);
    h_l = (ke_l + ue_l + p_l) / rho_l;

    rho_r = right_state[0];
    vel_r = right_state[1];
    p_r = right_state[3];
    eta_r = right_state[4];

    ke_r = 0.5 * rho_r * vel_r * vel_r;
    ue_r = p_r / (gamma - 1.0);
    h_r = (ke_r + ue_r + p_r) / rho_r;

    // roe averages
    Set::Scalar rho_RL, vel_RL, h_RL, a_RL;
    rho_RL = std::sqrt(rho_l) * std::sqrt(rho_r);
    vel_RL = (std::sqrt(rho_r) * vel_r + std::sqrt(rho_l) * vel_l) / (std::sqrt(rho_l) + std::sqrt(rho_r));
    h_RL = (std::sqrt(rho_r) * h_r + std::sqrt(rho_l) * h_l) / (std::sqrt(rho_l) + std::sqrt(rho_r));
    a_RL = std::sqrt((gamma - 1) * (h_RL - 0.5 * vel_RL * vel_RL));

    // eigenvalues
    std::array<Set::Scalar, 4> ev;
    ev[0] = vel_RL;
    ev[1] = vel_RL + a_RL;
    ev[2] = 0;
    ev[3] = vel_RL - a_RL;

    std::array<Set::Scalar, 4> dv;
    dv[0] = (rho_r - rho_l) - (p_r - p_l) / (a_RL * a_RL);
    dv[1] = (vel_r - vel_l) + (p_r - p_l) / (rho_RL * a_RL);
    dv[2] = 0;
    dv[2] = (vel_r - vel_l) - (p_r - p_l) / (rho_RL * a_RL);

    // right_state characteristic eigenvectors
    //std::array<Set::Scalar,4> rev;
    Eigen::Matrix<Set::Scalar, 4, 4> rev;
    rev(0, 0) = 1;
    rev(0, 1) = vel_RL;
    rev(0, 2) = 0;
    rev(0, 3) = 0.5 * vel_RL * vel_RL;

    rev(1, 0) = (rho_RL / 2 / a_RL);
    rev(1, 1) = (rho_RL / 2 / a_RL) * (vel_RL + a_RL);
    rev(1, 2) = 0;
    rev(1, 3) = (rho_RL / 2 / a_RL) * (h_RL + a_RL * vel_RL);

    rev(1, 0) = 0;
    rev(1, 1) = 0;
    rev(1, 2) = 0;
    rev(1, 3) = 0;

    rev(2, 0) = (-rho_RL / 2 / a_RL);
    rev(2, 1) = (-rho_RL / 2 / a_RL) * (vel_RL - a_RL);
    rev(2, 2) = 0;
    rev(2, 3) = (-rho_RL / 2 / a_RL) * (h_RL - a_RL * vel_RL);

    // compute fluxes
    std::array<Set::Scalar, 4> dynamic_flux_left;
    dynamic_flux_left[0] = left_state[0] * left_state[1];
    dynamic_flux_left[2] = dynamic_flux_left[0] * left_state[1];
    dynamic_flux_left[3] = dynamic_flux_left[0] * left_state[1];
    dynamic_flux_left[1] = left_state[1] * ke_l;

    std::array<Set::Scalar, 4> static_flux_left;
    static_flux_left[0] = 0.;
    static_flux_left[2] = left_state[3];
    static_flux_left[3] = 0.;
    static_flux_left[1] = left_state[1] * (ue_l + left_state[3]);

    std::array<Set::Scalar, 4> dynamic_flux_right;
    dynamic_flux_right[0] = right_state[0] * right_state[1];
    dynamic_flux_right[2] = dynamic_flux_right[0] * right_state[1];
    dynamic_flux_right[3] = dynamic_flux_right[0] * right_state[1];
    dynamic_flux_right[1] = right_state[1] * ke_r;

    std::array<Set::Scalar, 4> static_flux_right;
    static_flux_right[0] = 0.;
    static_flux_right[2] = right_state[3];
    static_flux_right[3] = 0.;
    static_flux_right[1] = right_state[1] * (ue_r + right_state[3]);

    // flux
    std::array<Set::Scalar, 4> fl;
    std::array<Set::Scalar, 4> upwind;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            fl[j] = 0.5 * (dynamic_flux_left[j] * eta_l * left_state[4] + dynamic_flux_right[j] * eta_r * right_state[4] + 0.5 * (static_flux_left[j] + static_flux_right[j]) * eta);
            upwind[j] = -0.5 * rev(i, j) * std::abs(ev[i]) * dv[i];
            fl[j] += upwind[j];
        };

    };

    return fl;
}
}
}
}
