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
// Every term vanishes when either phase is absent, so removing a phase from a
// calculation leaves the lower-dimensional model unchanged.  On a binary
// interface the equilibrium profile is
//
//   eta = 0.5 [1 - tanh(2 x / ell)]
//
// and the integrated interfacial energy is exactly sigma_ij.  This permits an
// independent surface tension for every material pair rather than restricting
// an N-phase system to additive per-phase surface energies.
//
// Liquid identities remain conserved partial-density fields.  This model
// supplies their thermodynamic free energy and chemical potentials only;
// actual material phase change remains the responsibility of a configured
// thermochemical mechanism.
//

#ifndef MODEL_PHASEFIELD_MULTLIQUID_H
#define MODEL_PHASEFIELD_MULTLIQUID_H

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
                const Set::Scalar interface_thickness,
                const std::vector<Set::Scalar>& surface_tension)
    {
        Util::Assert(INFO, TEST(names.size() >= 2));
        Util::Assert(INFO, TEST(surface_tension.size() ==
                                names.size() * names.size()));
        phase_names = names;
        epsilon = interface_thickness;
        sigma = surface_tension;
        sigma_max = 0.0;
        for (int i = 0; i < NumberOfPhases(); ++i)
            for (int j = i + 1; j < NumberOfPhases(); ++j)
                sigma_max = Util::Max(sigma_max, SurfaceTension(i,j));
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
    Set::Scalar MaximumSurfaceTension() const { return sigma_max; }

    Set::Scalar SurfaceTension(const int a, const int b) const
    {
        return sigma[a * NumberOfPhases() + b];
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

private:
    std::vector<std::string> phase_names;
    std::vector<Set::Scalar> sigma;
    Set::Scalar epsilon = 0.0;
    Set::Scalar sigma_max = 0.0;
};
}
}

#endif
