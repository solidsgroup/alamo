//
// Reduction-consistent diffuse-interface model for an arbitrary number of
// immiscible phases.  For phase fractions eta_i on the Gibbs simplex,
//
//   sum_i eta_i = 1,  eta_i >= 0,
//
// the pairwise free energy is
//
//   F_ij = c_ij [
//         3 ell/4 |eta_i grad(eta_j) - eta_j grad(eta_i)|^2
//       + 12/ell eta_i^2 eta_j^2 ].
//
// For liquid-liquid and liquid-gas pairs, c_ij is the physical surface
// tension.  A liquid-solid pair instead uses a positive coefficient kappa_ij
// in the same pair functional only to regularize the diffuse interface.  Its
// physical surface-energy difference relative to exposed solid is supplied by
//
//   F_surface = 2 (Delta_gamma_ij - kappa_ij)
//               h(eta_liquid) delta(eta_solid),
//
// where Delta_gamma_ij = gamma_SL - gamma_SG,
// h(q) = q^2(3-2q), and
//
//   delta(eta_solid) = sqrt(|grad(eta_solid)|^2 + r^2) - r.
//
// Since 2 integral_0^1 h(q) dq = 1, a flat binary liquid-solid
// interface has total energy kappa_ij + (Delta_gamma_ij-kappa_ij) =
// Delta_gamma_ij per unit area.  The subtraction of r makes the surface term
// exactly zero in uniform bulk regions.  h'(0) = h'(1) = 0 prevents a bulk
// force.
//
// A contact angle theta measured through the liquid is equivalent to
//
//   Delta_gamma_ij = -sigma_LG cos(theta)
//
// by Young's equation.  The input layer accepts either theta or Delta_gamma,
// never both, so they cannot independently over-constrain the same physics.
//
// On a binary interface the equilibrium profile is
//
//   eta = 0.5 [1 - tanh(2 x / ell)],
//
// and the integrated pair functional is exactly c_ij.  Every pair term is
// reduction-consistent: removing a phase leaves the lower-dimensional model
// unchanged.
//
// The phase list consists of every liquid, every solid, and an aggregate gas.
// Active interfaces must contain at least one liquid: liquid-liquid,
// liquid-gas, and optionally liquid-solid.  Solid-solid and solid-gas entries
// are always zero.  Thus Delta_gamma is a signed, referenced surface energy
// rather than an absolute solid surface tension.  Material identities remain
// conserved partial-density fields; actual phase change remains the
// responsibility of a thermochemical mechanism.
//

#ifndef MODEL_PHASEFIELD_MULTIPHASEINTERFACE_H
#define MODEL_PHASEFIELD_MULTIPHASEINTERFACE_H

#include <cmath>
#include <string>
#include <vector>

#include "AMReX_Array4.H"
#include "AMReX_Extension.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Model
{
namespace PhaseField
{
class MultiphaseInterface
{
public:
    static constexpr const char* name = "multiphase_interface";

    void Define(const std::vector<std::string>& names,
                const int number_liquid_phases,
                const bool include_solid_interfaces,
                const Set::Scalar interface_thickness,
                const std::vector<Set::Scalar>& surface_tension,
                const std::vector<Set::Scalar>&
                    solid_interface_regularization,
                const std::vector<Set::Scalar>&
                    solid_surface_energy_difference)
    {
        Util::Assert(INFO, TEST(names.size() >= 2));
        Util::Assert(INFO, TEST(number_liquid_phases >= 1));
        Util::Assert(INFO, TEST(number_liquid_phases <
                                static_cast<int>(names.size())));
        Util::Assert(INFO, TEST(interface_thickness > 0.0));
        Util::Assert(INFO, TEST(surface_tension.size() ==
                                names.size() * names.size()));
        Util::Assert(INFO, TEST(solid_interface_regularization.size() ==
                                names.size() * names.size()));
        Util::Assert(INFO, TEST(solid_surface_energy_difference.size() ==
                                names.size() * names.size()));
        phase_names = names;
        nliquid = number_liquid_phases;
        epsilon = interface_thickness;
        sigma = surface_tension;
        solid_regularization = solid_interface_regularization;
        solid_surface_difference = solid_surface_energy_difference;
        capillary_energy_max = 0.0;
        for (int i = 0; i < NumberOfPhases(); ++i)
            for (int j = i + 1; j < NumberOfPhases(); ++j)
            {
                const bool one_liquid =
                    IsLiquidPhase(i) || IsLiquidPhase(j);
                const bool one_solid = IsSolidPhase(i) || IsSolidPhase(j);
                const bool liquid_solid = one_liquid && one_solid;
                const bool fluid_interface = one_liquid && !one_solid;
                Util::Assert(INFO, TEST(fluid_interface ?
                    SurfaceTension(i,j) > 0.0 :
                    SurfaceTension(i,j) == 0.0));
                Util::Assert(INFO, TEST(SurfaceTension(i,j) ==
                                        SurfaceTension(j,i)));

                const bool active_solid = liquid_solid &&
                    include_solid_interfaces;
                Util::Assert(INFO, TEST(active_solid ?
                    SolidInterfaceRegularization(i,j) > 0.0 :
                    SolidInterfaceRegularization(i,j) == 0.0));
                Util::Assert(INFO, TEST(SolidInterfaceRegularization(i,j) ==
                    SolidInterfaceRegularization(j,i)));
                Util::Assert(INFO, TEST(active_solid ?
                    std::isfinite(SolidSurfaceEnergyDifference(i,j)) :
                    SolidSurfaceEnergyDifference(i,j) == 0.0));
                Util::Assert(INFO, TEST(SolidSurfaceEnergyDifference(i,j) ==
                    SolidSurfaceEnergyDifference(j,i)));

                Set::Scalar stiffness = SurfaceTension(i,j);
                if (active_solid)
                    stiffness = SolidInterfaceRegularization(i,j) +
                        2.0 * std::abs(SolidSurfaceCorrection(i,j));
                capillary_energy_max = Util::Max(
                    capillary_energy_max, stiffness);
            }
    }

    int NumberOfPhases() const
    {
        return static_cast<int>(phase_names.size());
    }

    const std::string& PhaseName(const int phase) const
    {
        return phase_names[phase];
    }

    Set::Scalar InterfaceThickness() const { return epsilon; }
    Set::Scalar MaximumCapillaryEnergy() const
    {
        return capillary_energy_max;
    }
    bool IsLiquidPhase(const int phase) const { return phase < nliquid; }
    bool IsSolidPhase(const int phase) const
    {
        return phase >= nliquid && phase < NumberOfPhases() - 1;
    }

    bool HasInterfacialEnergy(const int a, const int b) const
    {
        return PairRegularization(a,b) > 0.0;
    }

    Set::Scalar SurfaceTension(const int a, const int b) const
    {
        return sigma[a * NumberOfPhases() + b];
    }

    Set::Scalar SolidInterfaceRegularization(
        const int a, const int b) const
    {
        return solid_regularization[a * NumberOfPhases() + b];
    }

    Set::Scalar SolidSurfaceEnergyDifference(
        const int a, const int b) const
    {
        return solid_surface_difference[a * NumberOfPhases() + b];
    }

    Set::Scalar SolidSurfaceCorrection(const int a, const int b) const
    {
        return SolidSurfaceEnergyDifference(a,b) -
            SolidInterfaceRegularization(a,b);
    }

    Set::Scalar PairRegularization(const int a, const int b) const
    {
        return SurfaceTension(a,b) + SolidInterfaceRegularization(a,b);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar Interpolation(const Set::Scalar q)
    {
        const Set::Scalar bounded = Util::Clamp(q, 0.0, 1.0);
        return bounded * bounded * (3.0 - 2.0 * bounded);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar InterpolationDerivative(const Set::Scalar q)
    {
        const Set::Scalar bounded = Util::Clamp(q, 0.0, 1.0);
        return 6.0 * bounded * (1.0 - bounded);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar RegularizedSurfaceDelta(
        const Set::Scalar gradient_norm,
        const Set::Scalar interface_thickness)
    {
        const Set::Scalar regularization =
            1.0e-12 / interface_thickness;
        return std::sqrt(gradient_norm * gradient_norm +
                         regularization * regularization) - regularization;
    }

    static Set::Scalar SurfaceEnergyDifferenceFromContactAngle(
        const Set::Scalar liquid_gas_surface_tension,
        const Set::Scalar contact_angle)
    {
        return -liquid_gas_surface_tension * std::cos(contact_angle);
    }

    // Symmetric, pressure-reduced Korteweg stress associated with the pair
    // free energy.  The omitted f I term is an isotropic pressure gauge; its
    // discrete cancellation is handled by the pressure projection rather than
    // represented as a second, separately differenced pressure.  Applying
    // capillarity as the conservative divergence of this
    // stress makes every internal face traction equal and opposite and, since
    // the stress is symmetric, introduces no artificial capillary torque.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Matrix PairCapillaryStress(
        const amrex::Array4<const Set::Scalar>& eta,
        const int phase, const int other,
        const Set::Scalar pair_regularization,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar eta_i = eta(i,j,k,phase);
        const Set::Scalar eta_j = eta(i,j,k,other);
        const Set::Vector grad_i = Numeric::Gradient(
            eta, i, j, k, phase, dx, stencil);
        const Set::Vector grad_j = Numeric::Gradient(
            eta, i, j, k, other, dx, stencil);
        const Set::Vector relative_gradient =
            eta_i * grad_j - eta_j * grad_i;
        return -1.5 * interface_thickness * pair_regularization *
            relative_gradient * relative_gradient.transpose();
    }

    // Korteweg stress of the signed liquid-covered-solid surface correction
    //
    //   2 Delta h(eta_l) [sqrt(|grad eta_s|^2+r^2)-r].
    //
    // This is kept separate from PairCapillaryStress because the positive
    // pair coefficient regularizes the liquid-solid profile whereas Delta
    // selects the physical solid surface-energy difference/contact angle.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Matrix SolidSurfaceCorrectionCapillaryStress(
        const amrex::Array4<const Set::Scalar>& eta,
        const int liquid, const int solid,
        const Set::Scalar surface_correction,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar h = Interpolation(eta(i,j,k,liquid));
        const Set::Vector grad_solid = Numeric::Gradient(
            eta, i, j, k, solid, dx, stencil);
        const Set::Scalar regularization =
            1.0e-12 / interface_thickness;
        const Set::Scalar regularized_norm = std::sqrt(
            grad_solid.dot(grad_solid) +
            regularization * regularization);
        return -2.0 * surface_correction * h / regularized_norm *
            grad_solid * grad_solid.transpose();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar PairChemicalPotential(
        const amrex::Array4<const Set::Scalar>& eta,
        const int phase, const int other,
        const Set::Scalar pair_regularization,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar eta_i = eta(i,j,k,phase);
        const Set::Scalar eta_j = eta(i,j,k,other);
        const Set::Vector grad_i = Numeric::Gradient(
            eta, i, j, k, phase, dx, stencil);
        const Set::Vector grad_j = Numeric::Gradient(
            eta, i, j, k, other, dx, stencil);
        const Set::Scalar lap_i = Numeric::Laplacian(
            eta, i, j, k, phase, dx, stencil);
        const Set::Scalar lap_j = Numeric::Laplacian(
            eta, i, j, k, other, dx, stencil);
        const Set::Vector relative_gradient =
            eta_i * grad_j - eta_j * grad_i;
        const Set::Scalar relative_laplacian =
            eta_i * lap_j - eta_j * lap_i;

        return 3.0 * interface_thickness * pair_regularization *
                    relative_gradient.dot(grad_j) +
               1.5 * interface_thickness * pair_regularization * eta_j *
                    relative_laplacian +
               24.0 * pair_regularization / interface_thickness *
                    eta_i * eta_j * eta_j;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar LiquidSurfaceCorrectionChemicalPotential(
        const amrex::Array4<const Set::Scalar>& eta,
        const int liquid, const int solid,
        const Set::Scalar surface_correction,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar hprime = InterpolationDerivative(
            eta(i,j,k,liquid));
        const Set::Vector grad_solid = Numeric::Gradient(
            eta, i, j, k, solid, dx, stencil);
        const Set::Scalar surface_delta = RegularizedSurfaceDelta(
            std::sqrt(grad_solid.dot(grad_solid)), interface_thickness);
        return 2.0 * surface_correction * hprime * surface_delta;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar SolidSurfaceCorrectionChemicalPotential(
        const amrex::Array4<const Set::Scalar>& eta,
        const int liquid, const int solid,
        const Set::Scalar surface_correction,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar q = eta(i,j,k,liquid);
        const Set::Scalar h = Interpolation(q);
        const Set::Scalar hprime = InterpolationDerivative(q);
        const Set::Vector grad_liquid = Numeric::Gradient(
            eta, i, j, k, liquid, dx, stencil);
        const Set::Vector grad_solid = Numeric::Gradient(
            eta, i, j, k, solid, dx, stencil);
        const Set::Matrix hessian_solid = Numeric::Hessian(
            eta, i, j, k, solid, dx, stencil);
        const Set::Scalar regularization = 1.0e-12 / interface_thickness;
        const Set::Scalar surface_delta = std::sqrt(
            grad_solid.dot(grad_solid) +
            regularization * regularization);
        const Set::Vector normal = grad_solid / surface_delta;
        Set::Scalar projected_hessian = 0.0;
        Set::Scalar laplacian = 0.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            laplacian += hessian_solid(d,d);
            for (int e = 0; e < AMREX_SPACEDIM; ++e)
                projected_hessian += grad_solid(d) *
                    hessian_solid(d,e) * grad_solid(e);
        }
        const Set::Scalar curvature = laplacian / surface_delta -
            projected_hessian /
                (surface_delta * surface_delta * surface_delta);
        return -2.0 * surface_correction *
            (hprime * grad_liquid.dot(normal) + h * curvature);
    }

private:
    std::vector<std::string> phase_names;
    std::vector<Set::Scalar> sigma;
    std::vector<Set::Scalar> solid_regularization;
    std::vector<Set::Scalar> solid_surface_difference;
    int nliquid = 0;
    Set::Scalar epsilon = 0.0;
    Set::Scalar capillary_energy_max = 0.0;
};
}
}

#endif
