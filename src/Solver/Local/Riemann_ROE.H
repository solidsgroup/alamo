#include <eigen3/Eigen/Core>
#include "Set/Set.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{
/// Roe Riemann Solver based on Gas Dynamics - Culbert B. Laney
  void Riemann_ROE(double left[5], double right[5], eta) {
    Set::Scalar rho_l, vel_l, p_l, ke_l, ue_l, h_l;
    Set::Scalar rho_r, vel_r, p_r, ke_r, ue_r, h_r;
    Set::Scalar gamma=1.4;
    
    rho_l = left[0];
    vel_l = left[1];
    p_l = left[3];
    ke_l = 0.5 * rho_l * vel_l * vel_l;
    ue_l = p_l / (gamma - 1.0);
    h_l = (ke_l + ue_l + p_l) / rho_l;
        
    rho_r = right[0];
    vel_r = right[1];
    p_r = right[3];
    ke_r = 0.5 * rho_r * vel_r * vel_r;
    ue_r = p_r / (gamma - 1.0);
    h_r = (ke_r + ue_r + p_r) / rho_r;
        
    // roe averages
    Set::Scalar rho_RL, vel_RL, h_RL, a_RL;
    rho_RL = std::sqrt(rho_l) * std::sqrt(rho_r);
    vel_RL = (std::sqrt(rho_r) * vel_r + std::sqrt(rho_l) * vel_l) / (std::sqrt(rho_l) + std::sqrt(rho_r));
    h_RL = (std::sqrt(rho_r) * h_r + std::sqrt(rho_l) * h_l) / (std::sqrt(rho_l) + std::sqrt(rho_r));
    a_RL = std::sqrt((gamma - 1) * (h_RL - 0.5 * vel_RL ** 2));
        
    // eigenvalues
    Set::Matrix ev = Set::Matrix::Zero();
    ev[0] = vel_RL;
    ev[1] = vel_RL + a_RL;
    ev[2] = vel_RL - a_RL;

    Set::Matrix dv = Set::Matrix::Zero();
    dv[0] = (rho_r - rho_l) - (p_r - p_l) / a_RL2;
    dv[1] = (vel_r - vel_l) + (p_r - p_l) / (rho_RL * a_RL);
    dv[2] = (vel_r - vel_l) - (p_r - p_l) / (rho_RL * a_RL);
        
    // right characteristic eigenvectors
    Set::Matrix rev = Set::Matrix::Zero();
    rev[0, 0] = 1;
    rev[0, 1] = vel_RL;
    rev[0, 2] = 0.5 * vel_RL ** 2;

    rev[1, 0] = (rho_RL / 2 / a_RL);
    rev[1, 1] = (rho_RL / 2 / a_RL) * (vel_RL + a_RL);
    rev[1, 2] = (rho_RL / 2 / a_RL) * (h_RL + a_RL * vel_RL);

    rev[2, 0] = (-rho_RL / 2 / a_RL);
    rev[2, 1] = (-rho_RL / 2 / a_RL) * (vel_RL - a_RL);
    rev[2, 2] = (-rho_RL / 2 / a_RL) * (h_RL - a_RL * vel_RL);

    // compute fluxes
    Set::Scalar ue_left, ke_left, ue_right, ke_right;
    ue_left = ue_l;
    ke_left = 0.5 * left[0] * (left[1] * left[1] + left[2] * left[2]);

    ue_right = ue_r;
    ke_right = 0.5 * left[0] * (left[1] * left[1] + left[2] * left[2]);

    Set::Matrix dynamic_flux_left = Set::Matrix::Zero(); 
    dynamic_flux_left[0] = left[0] * left[1];
    dynamic_flux_left[2] = dynamic_flux_left[0] * left[1];
    dynamic_flux_left[3] = dynamic_flux_left[0] * left[1];
    dynamic_flux_left[1] = left[1] * ke_left;

    Set::Matrix static_flux_left = Set::Matrix::Zero(); 
    static_flux_left[0] = 0.;
    static_flux_left[2] = left[3];
    static_flux_left[3] = 0.;
    static_flux_left[1] = left[1] * (ue_left + left[3]);

    Set::Matrix dynamic_flux_right = Set::Matrix::Zero(); 
    dynamic_flux_right[0] = right[0] * right[1];
    dynamic_flux_right[2] = dynamic_flux_right[0] * right[1];
    dynamic_flux_right[3] = dynamic_flux_right[0] * right[1];
    dynamic_flux_right[1] = right[1] * ke_right;

    Set::Matrix static_flux_right = Set::Matrix::Zero(); 
    static_flux_right[0] = 0.;
    static_flux_right[2] = right[3];
    static_flux_right[3] = 0.;
    static_flux_right[1] = right[1] * (ue_right + right[3]);
        
    // flux
    Set::Matrix fl = Set::Matrix::Zero();
    Set::Matrix upwind = Set::Matrix::Zero();
    for (int i = 0; i < 4; i++) {
      fl[i] = 0.5 * (dynamic_flux_left[i] * left[4] + dynamic_flux_right[i] * right[4] + 0.5 * (static_flux_left[i] + static_flux_right[i]) * eta);
      for (int j = 0; j < 4; j++) {
	upwind[i,j] = -0.5 * rev[i, j] * std::abs(ev[i]) * dv[i];
	fl[i] = fl[i] + upwind[i,j]
      };
    
    };

    return fl;
}
