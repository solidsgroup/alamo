//
// Reduction-consistent diffuse-interface model for an arbitrary number of
// immiscible phases.  For phase fractions eta_i on the Gibbs simplex,
//
//   sum_i eta_i = 1,  eta_i >= 0,
//
// the pairwise free energy is
//
//   F = sum_{i<j} sigma_ij [
//         3 ell/4 |eta_i grad(eta_j) - eta_j grad(eta_i)|^2
//       + 12/ell eta_i^2 eta_j^2 ].
//
// A liquid-solid pair may additionally have the solid-surface energy
//
//       -A_ij h(eta_liquid) |grad(eta_solid)|,
//
// where h(q) = q^2(3-2q).  It lowers the energy by A_ij per unit solid surface
// covered by liquid, while h'(0) = h'(1) = 0 prevents a bulk force.  Unlike a
// bilinear eta_liquid*eta_solid term, its variational force is not a pure
// pressure gradient.  It also assigns no energy to an exposed solid-gas
// interface.
//
// Every term vanishes when either phase is absent, so removing a phase from a
// calculation leaves the lower-dimensional model unchanged.  On a binary
// interface the equilibrium profile is
//
//   eta = 0.5 [1 - tanh(2 x / ell)]
//
// and the integrated interfacial energy is exactly sigma_ij.  This permits an
// independent interfacial energy for every active material pair rather than
// restricting an N-phase system to additive per-phase surface energies.
//
// The phase list consists of every liquid, every solid, and an aggregate gas.
// Active interfaces must contain at least one liquid: liquid-liquid,
// liquid-gas, and optionally liquid-solid.  Solid-solid and solid-gas entries
// are always zero, so adding liquid capillarity does not introduce an
// independent solid surface model.  Material identities remain conserved
// partial-density fields; actual phase change remains the responsibility of a
// thermochemical mechanism.
//

#ifndef MODEL_PHASEFIELD_MULTLIQUID_H
#define MODEL_PHASEFIELD_MULTLIQUID_H

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
class MultLiquid
{
public:
    static constexpr const char* name = "multliquid";

    void Define(const std::vector<std::string>& names,
                const int number_liquid_phases,
                const bool include_solid_interfaces,
                const Set::Scalar interface_thickness,
                const std::vector<Set::Scalar>& surface_tension,
                const std::vector<Set::Scalar>& adhesion_energy)
    {
        Util::Assert(INFO, TEST(names.size() >= 2));
        Util::Assert(INFO, TEST(number_liquid_phases >= 1));
        Util::Assert(INFO, TEST(number_liquid_phases <
                                static_cast<int>(names.size())));
        Util::Assert(INFO, TEST(surface_tension.size() ==
                                names.size() * names.size()));
        Util::Assert(INFO, TEST(adhesion_energy.size() ==
                                names.size() * names.size()));
        phase_names = names;
        nliquid = number_liquid_phases;
        epsilon = interface_thickness;
        sigma = surface_tension;
        adhesion = adhesion_energy;
        capillary_energy_max = 0.0;
        for (int i = 0; i < NumberOfPhases(); ++i)
            for (int j = i + 1; j < NumberOfPhases(); ++j)
            {
                const bool one_liquid =
                    IsLiquidPhase(i) || IsLiquidPhase(j);
                const bool one_solid = IsSolidPhase(i) || IsSolidPhase(j);
                const bool active = one_liquid &&
                    (include_solid_interfaces || !one_solid);
                Util::Assert(INFO, TEST(active ? SurfaceTension(i,j) > 0.0 :
                                                 SurfaceTension(i,j) == 0.0));
                Util::Assert(INFO, TEST(SurfaceTension(i,j) ==
                                        SurfaceTension(j,i)));
                const bool liquid_solid = one_liquid && one_solid;
                Util::Assert(INFO, TEST(liquid_solid ?
                    AdhesionEnergy(i,j) >= 0.0 :
                    AdhesionEnergy(i,j) == 0.0));
                Util::Assert(INFO, TEST(AdhesionEnergy(i,j) ==
                                        AdhesionEnergy(j,i)));
                capillary_energy_max = Util::Max(capillary_energy_max,
                    SurfaceTension(i,j) + AdhesionEnergy(i,j));
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
        return SurfaceTension(a,b) > 0.0;
    }

    Set::Scalar SurfaceTension(const int a, const int b) const
    {
        return sigma[a * NumberOfPhases() + b];
    }

    Set::Scalar AdhesionEnergy(const int a, const int b) const
    {
        return adhesion[a * NumberOfPhases() + b];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar PairChemicalPotential(
        const amrex::Array4<const Set::Scalar>& eta,
        const int phase, const int other,
        const Set::Scalar surface_tension,
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

        return 3.0 * interface_thickness * surface_tension *
                    relative_gradient.dot(grad_j) +
               1.5 * interface_thickness * surface_tension * eta_j *
                    relative_laplacian +
               24.0 * surface_tension / interface_thickness *
                    eta_i * eta_j * eta_j;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar LiquidAdhesionChemicalPotential(
        const amrex::Array4<const Set::Scalar>& eta,
        const int liquid, const int solid,
        const Set::Scalar adhesion_energy,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar q = Util::Clamp(eta(i,j,k,liquid), 0.0, 1.0);
        const Set::Scalar hprime = 6.0 * q * (1.0 - q);
        const Set::Vector grad_solid = Numeric::Gradient(
            eta, i, j, k, solid, dx, stencil);
        const Set::Scalar regularization = 1.0e-12 / interface_thickness;
        const Set::Scalar surface_delta = std::sqrt(
            grad_solid.dot(grad_solid) +
            regularization * regularization);
        return -adhesion_energy * hprime * surface_delta;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar SolidAdhesionChemicalPotential(
        const amrex::Array4<const Set::Scalar>& eta,
        const int liquid, const int solid,
        const Set::Scalar adhesion_energy,
        const Set::Scalar interface_thickness,
        const int i, const int j, const int k,
        const Set::Scalar dx[AMREX_SPACEDIM],
        std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil)
    {
        const Set::Scalar q = Util::Clamp(eta(i,j,k,liquid), 0.0, 1.0);
        const Set::Scalar h = q * q * (3.0 - 2.0 * q);
        const Set::Scalar hprime = 6.0 * q * (1.0 - q);
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
        return adhesion_energy *
            (hprime * grad_liquid.dot(normal) + h * curvature);
    }

private:
    std::vector<std::string> phase_names;
    std::vector<Set::Scalar> sigma;
    std::vector<Set::Scalar> adhesion;
    int nliquid = 0;
    Set::Scalar epsilon = 0.0;
    Set::Scalar capillary_energy_max = 0.0;
};
}
}

#endif
