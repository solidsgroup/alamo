#include "LowMach.H"

#include "AMReX_MLABecLaplacian.H"
#include "AMReX_FillPatchUtil.H"
#include "AMReX_MLMG.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_TimeIntegrator.H"
#include "Numeric/Stencil.H"

#include "Model/Gas/Thermo/Thermo.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Transport/Transport.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/EOS/EOS.H"
#include "Model/Gas/EOS/CPG.H"

namespace
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldWeight(const Set::Scalar eta_val, const Set::Scalar eta_band)
{
    if (eta_val <= eta_band || eta_val >= 1.0 - eta_band) return 0.0;
    return eta_val * (1.0 - eta_val);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldCoordinate(const Set::Scalar eta_val,
                        const Set::Scalar epsilon,
                        const Set::Scalar eta_band,
                        const bool mapped_distance)
{
    if (!mapped_distance) return eta_val;

    // Recover the distance coordinate of the target tanh profile.
    const Set::Scalar eta_floor = Util::Clamp(0.1 * eta_band, 1.0e-12, 0.25);
    const Set::Scalar bounded_eta = Util::Clamp(eta_val, eta_floor, 1.0 - eta_floor);
    return 0.5 * epsilon * std::log(bounded_eta / (1.0 - bounded_eta));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Vector
EtaPhaseFieldCellNormal(const amrex::Array4<const Set::Scalar>& eta,
                        const int i, const int j, const int k,
                        const Set::Scalar DX[AMREX_SPACEDIM],
                        const Set::Scalar epsilon,
                        const Set::Scalar eta_band,
                        const Set::Scalar normal_regularization)
{
    Set::Vector gradient = Set::Vector::Zero();
    gradient(0) =
        (EtaPhaseFieldCoordinate(eta(i+1,j,k), epsilon, eta_band, true) -
         EtaPhaseFieldCoordinate(eta(i-1,j,k), epsilon, eta_band, true)) / (2.0 * DX[0]);
#if AMREX_SPACEDIM > 1
    gradient(1) =
        (EtaPhaseFieldCoordinate(eta(i,j+1,k), epsilon, eta_band, true) -
         EtaPhaseFieldCoordinate(eta(i,j-1,k), epsilon, eta_band, true)) / (2.0 * DX[1]);
#endif
#if AMREX_SPACEDIM > 2
    gradient(2) =
        (EtaPhaseFieldCoordinate(eta(i,j,k+1), epsilon, eta_band, true) -
         EtaPhaseFieldCoordinate(eta(i,j,k-1), epsilon, eta_band, true)) / (2.0 * DX[2]);
#endif
    const Set::Scalar norm =
        std::sqrt(gradient.squaredNorm() + normal_regularization * normal_regularization);
    gradient /= norm;
    return gradient;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldNormalCoherence(const amrex::Array4<const Set::Scalar>& eta,
                             const int i, const int j, const int k,
                             const int face_direction,
                             const Set::Scalar DX[AMREX_SPACEDIM],
                             const Set::Scalar epsilon,
                             const Set::Scalar eta_band,
                             const Set::Scalar normal_regularization,
                             const Set::Scalar threshold,
                             const Set::Scalar power)
{
    if (threshold <= 0.0) return 1.0;
#if AMREX_SPACEDIM == 1
    return 1.0;
#else

    Set::Vector normal_sum = Set::Vector::Zero();
    int count = 0;
    for (int side = 0; side <= 1; ++side)
    {
#if AMREX_SPACEDIM == 2
        const int tangent = 1 - face_direction;
        for (int offset = -1; offset <= 1; ++offset)
        {
            int index[3] = {i, j, k};
            index[face_direction] += side;
            index[tangent] += offset;
            normal_sum += EtaPhaseFieldCellNormal(eta, index[0], index[1], index[2], DX,
                                                   epsilon, eta_band, normal_regularization);
            ++count;
        }
#else
        const int tangent0 = (face_direction + 1) % 3;
        const int tangent1 = (face_direction + 2) % 3;
        for (int offset0 = -1; offset0 <= 1; ++offset0)
        for (int offset1 = -1; offset1 <= 1; ++offset1)
        {
            int index[3] = {i, j, k};
            index[face_direction] += side;
            index[tangent0] += offset0;
            index[tangent1] += offset1;
            normal_sum += EtaPhaseFieldCellNormal(eta, index[0], index[1], index[2], DX,
                                                   epsilon, eta_band, normal_regularization);
            ++count;
        }
#endif
    }

    const Set::Scalar coherence =
        std::sqrt(normal_sum.squaredNorm()) / static_cast<Set::Scalar>(count);
    // Grade reinitialization continuously where a sharp corner mixes normals.
    const Set::Scalar grade = Util::Clamp((coherence - threshold) / (1.0 - threshold), 0.0, 1.0);
    return std::pow(grade, power);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldCellCoherence(const amrex::Array4<const Set::Scalar>& eta,
                           const int i, const int j, const int k,
                           const Set::Scalar DX[AMREX_SPACEDIM],
                           const Set::Scalar epsilon,
                           const Set::Scalar eta_band,
                           const Set::Scalar normal_regularization,
                           const Set::Scalar threshold,
                           const Set::Scalar power)
{
    if (threshold <= 0.0) return 1.0;
#if AMREX_SPACEDIM == 1
    return 1.0;
#else
    Set::Vector normal_sum = Set::Vector::Zero();
    int count = 0;
#if AMREX_SPACEDIM == 2
    for (int ii = -1; ii <= 1; ++ii)
    for (int jj = -1; jj <= 1; ++jj)
    {
        normal_sum += EtaPhaseFieldCellNormal(eta, i+ii, j+jj, k, DX,
                                               epsilon, eta_band, normal_regularization);
        ++count;
    }
#else
    for (int ii = -1; ii <= 1; ++ii)
    for (int jj = -1; jj <= 1; ++jj)
    for (int kk = -1; kk <= 1; ++kk)
    {
        normal_sum += EtaPhaseFieldCellNormal(eta, i+ii, j+jj, k+kk, DX,
                                               epsilon, eta_band, normal_regularization);
        ++count;
    }
#endif
    const Set::Scalar coherence =
        std::sqrt(normal_sum.squaredNorm()) / static_cast<Set::Scalar>(count);
    const Set::Scalar grade = Util::Clamp((coherence - threshold) / (1.0 - threshold), 0.0, 1.0);
    return std::pow(grade, power);
#endif
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldFluxX(const amrex::Array4<const Set::Scalar>& eta,
                   const int i, const int j, const int k,
                   const Set::Scalar DX[AMREX_SPACEDIM],
                   const Set::Scalar epsilon,
                   const Set::Scalar counter_curvature,
                   const Set::Scalar eta_band,
                   const Set::Scalar normal_regularization,
                   const bool mapped_distance,
                   const Set::Scalar normal_coherence_threshold,
                   const Set::Scalar normal_coherence_power,
                   const Set::Scalar incoherent_diffusion)
{
    Set::Vector grad_coordinate = Set::Vector::Zero();
    grad_coordinate(0) =
        (EtaPhaseFieldCoordinate(eta(i+1,j,k), epsilon, eta_band, mapped_distance) -
         EtaPhaseFieldCoordinate(eta(i,j,k), epsilon, eta_band, mapped_distance)) / DX[0];
#if AMREX_SPACEDIM > 1
    grad_coordinate(1) = 0.25 *
        ((EtaPhaseFieldCoordinate(eta(i,j+1,k), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i,j-1,k), epsilon, eta_band, mapped_distance)) +
         (EtaPhaseFieldCoordinate(eta(i+1,j+1,k), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i+1,j-1,k), epsilon, eta_band, mapped_distance))) / DX[1];
#endif
#if AMREX_SPACEDIM > 2
    grad_coordinate(2) = 0.25 *
        ((EtaPhaseFieldCoordinate(eta(i,j,k+1), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i,j,k-1), epsilon, eta_band, mapped_distance)) +
         (EtaPhaseFieldCoordinate(eta(i+1,j,k+1), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i+1,j,k-1), epsilon, eta_band, mapped_distance))) / DX[2];
#endif
    const Set::Scalar eta_face = Util::Clamp(0.5 * (eta(i,j,k) + eta(i+1,j,k)), 0.0, 1.0);
    const Set::Scalar eta_weight = EtaPhaseFieldWeight(eta_face, eta_band);
    if (eta_weight == 0.0) return 0.0;

    const Set::Scalar grad_norm =
        std::sqrt(grad_coordinate.squaredNorm() + normal_regularization * normal_regularization);
    if (mapped_distance)
    {
        const Set::Scalar normal_coherence = EtaPhaseFieldNormalCoherence(
            eta, i, j, k, 0, DX, epsilon, eta_band, normal_regularization,
            normal_coherence_threshold, normal_coherence_power);
        const Set::Scalar diffusion_grade =
            normal_coherence + (1.0 - normal_coherence) * incoherent_diffusion;
        return eta_weight * grad_coordinate(0) *
               (diffusion_grade - normal_coherence * counter_curvature / grad_norm);
    }
    return 0.5 * epsilon * grad_coordinate(0)
           - counter_curvature * eta_weight * grad_coordinate(0) / grad_norm;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldFluxY(const amrex::Array4<const Set::Scalar>& eta,
                   const int i, const int j, const int k,
                   const Set::Scalar DX[AMREX_SPACEDIM],
                   const Set::Scalar epsilon,
                   const Set::Scalar counter_curvature,
                   const Set::Scalar eta_band,
                   const Set::Scalar normal_regularization,
                   const bool mapped_distance,
                   const Set::Scalar normal_coherence_threshold,
                   const Set::Scalar normal_coherence_power,
                   const Set::Scalar incoherent_diffusion)
{
    Set::Vector grad_coordinate = Set::Vector::Zero();
    grad_coordinate(0) = 0.25 *
        ((EtaPhaseFieldCoordinate(eta(i+1,j,k), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i-1,j,k), epsilon, eta_band, mapped_distance)) +
         (EtaPhaseFieldCoordinate(eta(i+1,j+1,k), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i-1,j+1,k), epsilon, eta_band, mapped_distance))) / DX[0];
#if AMREX_SPACEDIM > 1
    grad_coordinate(1) =
        (EtaPhaseFieldCoordinate(eta(i,j+1,k), epsilon, eta_band, mapped_distance) -
         EtaPhaseFieldCoordinate(eta(i,j,k), epsilon, eta_band, mapped_distance)) / DX[1];
#endif
#if AMREX_SPACEDIM > 2
    grad_coordinate(2) = 0.25 *
        ((EtaPhaseFieldCoordinate(eta(i,j,k+1), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i,j,k-1), epsilon, eta_band, mapped_distance)) +
         (EtaPhaseFieldCoordinate(eta(i,j+1,k+1), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i,j+1,k-1), epsilon, eta_band, mapped_distance))) / DX[2];
#endif
    const Set::Scalar eta_face = Util::Clamp(0.5 * (eta(i,j,k) + eta(i,j+1,k)), 0.0, 1.0);
    const Set::Scalar eta_weight = EtaPhaseFieldWeight(eta_face, eta_band);
    if (eta_weight == 0.0) return 0.0;

    const Set::Scalar grad_norm =
        std::sqrt(grad_coordinate.squaredNorm() + normal_regularization * normal_regularization);
    if (mapped_distance)
    {
        const Set::Scalar normal_coherence = EtaPhaseFieldNormalCoherence(
            eta, i, j, k, 1, DX, epsilon, eta_band, normal_regularization,
            normal_coherence_threshold, normal_coherence_power);
        const Set::Scalar diffusion_grade =
            normal_coherence + (1.0 - normal_coherence) * incoherent_diffusion;
        return eta_weight * grad_coordinate(1) *
               (diffusion_grade - normal_coherence * counter_curvature / grad_norm);
    }
    return 0.5 * epsilon * grad_coordinate(1)
           - counter_curvature * eta_weight * grad_coordinate(1) / grad_norm;
}

#if AMREX_SPACEDIM > 2
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldFluxZ(const amrex::Array4<const Set::Scalar>& eta,
                   const int i, const int j, const int k,
                   const Set::Scalar DX[AMREX_SPACEDIM],
                   const Set::Scalar epsilon,
                   const Set::Scalar counter_curvature,
                   const Set::Scalar eta_band,
                   const Set::Scalar normal_regularization,
                   const bool mapped_distance,
                   const Set::Scalar normal_coherence_threshold,
                   const Set::Scalar normal_coherence_power,
                   const Set::Scalar incoherent_diffusion)
{
    Set::Vector grad_coordinate = Set::Vector::Zero();
    grad_coordinate(0) = 0.25 *
        ((EtaPhaseFieldCoordinate(eta(i+1,j,k), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i-1,j,k), epsilon, eta_band, mapped_distance)) +
         (EtaPhaseFieldCoordinate(eta(i+1,j,k+1), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i-1,j,k+1), epsilon, eta_band, mapped_distance))) / DX[0];
    grad_coordinate(1) = 0.25 *
        ((EtaPhaseFieldCoordinate(eta(i,j+1,k), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i,j-1,k), epsilon, eta_band, mapped_distance)) +
         (EtaPhaseFieldCoordinate(eta(i,j+1,k+1), epsilon, eta_band, mapped_distance) -
          EtaPhaseFieldCoordinate(eta(i,j-1,k+1), epsilon, eta_band, mapped_distance))) / DX[1];
    grad_coordinate(2) =
        (EtaPhaseFieldCoordinate(eta(i,j,k+1), epsilon, eta_band, mapped_distance) -
         EtaPhaseFieldCoordinate(eta(i,j,k), epsilon, eta_band, mapped_distance)) / DX[2];
    const Set::Scalar eta_face = Util::Clamp(0.5 * (eta(i,j,k) + eta(i,j,k+1)), 0.0, 1.0);
    const Set::Scalar eta_weight = EtaPhaseFieldWeight(eta_face, eta_band);
    if (eta_weight == 0.0) return 0.0;

    const Set::Scalar grad_norm =
        std::sqrt(grad_coordinate.squaredNorm() + normal_regularization * normal_regularization);
    if (mapped_distance)
    {
        const Set::Scalar normal_coherence = EtaPhaseFieldNormalCoherence(
            eta, i, j, k, 2, DX, epsilon, eta_band, normal_regularization,
            normal_coherence_threshold, normal_coherence_power);
        const Set::Scalar diffusion_grade =
            normal_coherence + (1.0 - normal_coherence) * incoherent_diffusion;
        return eta_weight * grad_coordinate(2) *
               (diffusion_grade - normal_coherence * counter_curvature / grad_norm);
    }
    return 0.5 * epsilon * grad_coordinate(2)
           - counter_curvature * eta_weight * grad_coordinate(2) / grad_norm;
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
EtaPhaseFieldSource(const amrex::Array4<const Set::Scalar>& eta,
                    const int i, const int j, const int k,
                    const Set::Scalar DX[AMREX_SPACEDIM],
                    const Set::Scalar epsilon,
                    const Set::Scalar mobility,
                    const Set::Scalar counter_curvature,
                    const Set::Scalar eta_band,
                    const Set::Scalar normal_regularization,
                    const bool mapped_distance,
                    const Set::Scalar normal_coherence_threshold,
                    const Set::Scalar normal_coherence_power,
                    const Set::Scalar exterior_reinitialization,
                    const Set::Scalar incoherent_diffusion)
{
    if (mobility == 0.0 || epsilon <= 0.0) return 0.0;

    Set::Scalar div_flux =
        (EtaPhaseFieldFluxX(eta, i,   j, k, DX, epsilon, counter_curvature, eta_band, normal_regularization, mapped_distance,
                            normal_coherence_threshold, normal_coherence_power, incoherent_diffusion) -
         EtaPhaseFieldFluxX(eta, i-1, j, k, DX, epsilon, counter_curvature, eta_band, normal_regularization, mapped_distance,
                            normal_coherence_threshold, normal_coherence_power, incoherent_diffusion)) / DX[0];
#if AMREX_SPACEDIM > 1
    div_flux +=
        (EtaPhaseFieldFluxY(eta, i, j,   k, DX, epsilon, counter_curvature, eta_band, normal_regularization, mapped_distance,
                            normal_coherence_threshold, normal_coherence_power, incoherent_diffusion) -
         EtaPhaseFieldFluxY(eta, i, j-1, k, DX, epsilon, counter_curvature, eta_band, normal_regularization, mapped_distance,
                            normal_coherence_threshold, normal_coherence_power, incoherent_diffusion)) / DX[1];
#endif
#if AMREX_SPACEDIM > 2
    div_flux +=
        (EtaPhaseFieldFluxZ(eta, i, j, k,   DX, epsilon, counter_curvature, eta_band, normal_regularization, mapped_distance,
                            normal_coherence_threshold, normal_coherence_power, incoherent_diffusion) -
         EtaPhaseFieldFluxZ(eta, i, j, k-1, DX, epsilon, counter_curvature, eta_band, normal_regularization, mapped_distance,
                            normal_coherence_threshold, normal_coherence_power, incoherent_diffusion)) / DX[2];
#endif
    // Relax only the exterior distance profile; the eta=0.5 contour remains fixed.
    if (mapped_distance && exterior_reinitialization > 0.0)
    {
        const Set::Scalar eta_cell = Util::Clamp(eta(i,j,k), 0.0, 1.0);
        const Set::Scalar eta_weight = EtaPhaseFieldWeight(eta_cell, eta_band);
        if (eta_weight > 0.0 && eta_cell < 0.5)
        {
            const Set::Scalar coordinate =
                EtaPhaseFieldCoordinate(eta_cell, epsilon, eta_band, true);
            Set::Vector gradient = Set::Vector::Zero();
            gradient(0) =
                (EtaPhaseFieldCoordinate(eta(i+1,j,k), epsilon, eta_band, true) -
                 EtaPhaseFieldCoordinate(eta(i-1,j,k), epsilon, eta_band, true)) / (2.0 * DX[0]);
#if AMREX_SPACEDIM > 1
            gradient(1) =
                (EtaPhaseFieldCoordinate(eta(i,j+1,k), epsilon, eta_band, true) -
                 EtaPhaseFieldCoordinate(eta(i,j-1,k), epsilon, eta_band, true)) / (2.0 * DX[1]);
#endif
#if AMREX_SPACEDIM > 2
            gradient(2) =
                (EtaPhaseFieldCoordinate(eta(i,j,k+1), epsilon, eta_band, true) -
                 EtaPhaseFieldCoordinate(eta(i,j,k-1), epsilon, eta_band, true)) / (2.0 * DX[2]);
#endif
            const Set::Scalar gradient_norm = std::sqrt(
                gradient.squaredNorm() + normal_regularization * normal_regularization);
            const Set::Scalar sign_regularization = 0.5 * epsilon;
            const Set::Scalar smooth_sign = coordinate / std::sqrt(
                coordinate * coordinate + sign_regularization * sign_regularization);
            const Set::Scalar coherence = EtaPhaseFieldCellCoherence(
                eta, i, j, k, DX, epsilon, eta_band, normal_regularization,
                normal_coherence_threshold, normal_coherence_power);
            div_flux += exterior_reinitialization * (1.0 - coherence) *
                (2.0 * eta_weight / epsilon) * smooth_sign * (1.0 - gradient_norm);
        }
    }
    return mobility * div_flux;
}
}

namespace Integrator
{
LowMach::LowMach(IO::ParmParse& pp) : LowMach()
{
    pp_queryclass(*this);
}

void
LowMach::Parse(LowMach& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::LowMach::Parse");

    pp.query_required("cfl", value.cfl);
    pp.query_default("cfl_v", value.cfl_v, 1.0e100);
    pp.query_default("small", value.small, 1.0e-12);
    pp.query_default("density_floor", value.density_floor, value.small);
    pp.query_default("temperature_floor", value.temperature_floor, value.small);
    pp.query_default("pressure_floor", value.pressure_floor, value.small);
    pp.query_default("thermodynamic_pressure", value.thermodynamic_pressure, "100000.0_Pa", Unit::Pressure());
    pp.query_default("pressure_scale", value.pressure_scale, 1.0);
    pp.query_default("projection.enabled", value.projection_enabled, true);
    pp.query_default("projection.tol_rel", value.projection_tol_rel, 1.0e-11);
    pp.query_default("projection.tol_abs", value.projection_tol_abs, 1.0e-12);
    pp.query_default("projection.verbose", value.projection_verbose, 0);
    pp.query_default("projection.node_reconstruction_sweeps", value.projection_node_reconstruction_sweeps, 1);
    pp.query_default("projection.nodal_iterations", value.projection_nodal_iterations, 200);
    pp.query_default("diagnostics.interval", value.diagnostics_interval, 0);
    pp.query_default("diagnostics.extended_fields", value.diagnostics_extended_fields, false);
    pp.query_default("projection.update_pressure", value.projection_update_pressure, true);
    pp.query_default("projection.predictor_pressure_gradient", value.projection_predictor_pressure_gradient, false);
    pp.query_default("include_viscosity", value.include_viscosity, true);
    pp.query_default("include_conduction", value.include_conduction, true);
    pp.query_default("advect_temperature", value.advect_temperature, true);
    pp.query_default("eta.initial_value", value.eta_initial_value, 0.0);
    pp.query_default("eta.phase_field.enabled", value.eta_phase_field_enabled, false);
    pp.query_default("eta.phase_field.epsilon", value.eta_phase_field_epsilon, "0.0", Unit::Length());
    pp.query_default("eta.phase_field.mobility", value.eta_phase_field_mobility, "0.0", Unit::Velocity());
    pp.query_default("eta.phase_field.counter_curvature", value.eta_phase_field_counter_curvature, 1.0);
    pp.query_default("eta.phase_field.band", value.eta_phase_field_band, 1.0e-8);
    pp.query_default("eta.phase_field.normal_regularization", value.eta_phase_field_normal_regularization, 1.0e-12);
    pp.query_default("eta.phase_field.mapped_distance", value.eta_phase_field_mapped_distance, true);
    pp.query_default("eta.phase_field.normal_coherence_threshold", value.eta_phase_field_normal_coherence_threshold, 0.0);
    pp.query_default("eta.phase_field.normal_coherence_power", value.eta_phase_field_normal_coherence_power, 1.0);
    pp.query_default("eta.phase_field.exterior_reinitialization", value.eta_phase_field_exterior_reinitialization, 0.0);
    pp.query_default("eta.phase_field.incoherent_diffusion", value.eta_phase_field_incoherent_diffusion, 0.0);
    if (value.eta_phase_field_enabled && value.eta_phase_field_epsilon <= 0.0)
        Util::Exception(INFO, "eta.phase_field.epsilon must be positive when eta.phase_field.enabled=1");
    if (value.eta_phase_field_normal_coherence_threshold < 0.0 ||
        value.eta_phase_field_normal_coherence_threshold >= 1.0)
        Util::Exception(INFO, "eta.phase_field.normal_coherence_threshold must be in [0,1)");
    if (value.eta_phase_field_normal_coherence_power <= 0.0)
        Util::Exception(INFO, "eta.phase_field.normal_coherence_power must be positive");
    if (value.eta_phase_field_exterior_reinitialization < 0.0)
        Util::Exception(INFO, "eta.phase_field.exterior_reinitialization must be nonnegative");
    if (value.eta_phase_field_incoherent_diffusion < 0.0 ||
        value.eta_phase_field_incoherent_diffusion > 1.0)
        Util::Exception(INFO, "eta.phase_field.incoherent_diffusion must be in [0,1]");
    pp.query_default("reference_map.eta_cutoff", value.reference_map_eta_cutoff, 0.5);
    pp.query_default("reference_map.eta_core", value.reference_map_eta_core, value.reference_map_eta_cutoff);
    pp.query_default("reference_map.eta_extension", value.reference_map_eta_extension, 1.0e-3);
    pp.query_default("reference_map.extrapolation_sweeps", value.reference_map_extrapolation_sweeps, 4);
    pp.query_default("reference_map.reconstruction_alpha", value.reference_map_reconstruction_alpha, 1.0);
    pp.query_default("reference_map.reconstruction_power", value.reference_map_reconstruction_power, 1.0);
    pp.query_default("reference_map.affine_tolerance", value.reference_map_affine_tolerance, 1.0e-7);
    pp.query_default("reference_map.smoothing_sweeps", value.reference_map_smoothing_sweeps, 2);
    pp.query_default("reference_map.smoothing_alpha", value.reference_map_smoothing_alpha, 1.0);
    pp.query_default("reference_map.smoothing_power", value.reference_map_smoothing_power, 2.0);
    pp.query_default("reference_map.stress_smoothing_sweeps", value.reference_map_stress_smoothing_sweeps, 0);
    pp.query_default("reference_map.stress_smoothing_alpha", value.reference_map_stress_smoothing_alpha, 0.1);

    std::string solid_model_type;
    pp.query_default("solid.model.type", solid_model_type, "none");
    pp.query_default("solid.model.eta_threshold", value.finite_solid_eta_threshold, 0.5);
    pp.query_default("solid.model.J_floor", value.finite_solid_J_floor, 1.0e-6);
    pp.query_default("solid.model.viscosity", value.finite_solid_viscosity, "0.0", Unit::Pressure() * Unit::Time());
    pp.query_default("solid.model.bulk_viscosity", value.finite_solid_bulk_viscosity, "0.0", Unit::Pressure() * Unit::Time());
    pp.query_default("solid.model.interface_viscosity", value.finite_solid_interface_viscosity, "0.0", Unit::Pressure() * Unit::Time());
    pp.query_default("solid.model.stress_rhs_sign", value.finite_solid_deviatoric_stress_divergence_sign, 0.0);
    pp.query_default("solid.model.deviatoric_stress_divergence_sign",
                     value.finite_solid_deviatoric_stress_divergence_sign,
                     value.finite_solid_deviatoric_stress_divergence_sign);
    if (solid_model_type == "none")
    {
        value.finite_solid_enabled = false;
    }
    else if (solid_model_type == "finite.neohookean")
    {
        value.finite_solid_enabled = true;
        pp.queryclass<Model::Solid::Finite::NeoHookean>("solid.model.finite.neohookean", value.finite_solid_model);
    }
    else
    {
        Util::Exception(INFO, solid_model_type, " is not a valid LowMach solid.model.type");
    }

    pp.select<Numeric::Advect::MUSCL,
              Numeric::Advect::Upwind,
              Numeric::Advect::Centered,
              Numeric::Advect::QUICK,
              Numeric::Advect::WENO5>("advection",value.advect);
    if (value.advect.PhiLocation() != Set::HC::Cell ||
        value.advect.VelocityLocation() != Set::HC::Cell)
        Util::Exception(INFO, "LowMach currently requires cell-centered phi and velocity advection data");

    pp.query_default("velocity_refinement_criterion", value.velocity_refinement_criterion, 1.0e100);
    pp.query_default("pressure_refinement_criterion", value.pressure_refinement_criterion, 1.0e100);
    pp.query_default("temperature_refinement_criterion", value.temperature_refinement_criterion, 1.0e100);
    pp.query_default("eta_refinement_criterion", value.eta_refinement_criterion, 1.0e100);
    pp.queryarr_default("g", value.g, Set::Vector::Zero());

    pp.queryclass<Model::Gas::Gas>("gas", value.gas);
    value.nspecies = value.gas.nspecies;

    int nghost = value.advect.NGhost();
    if (nghost < 3) nghost = 3;
    pp.select_default<BC::Constant,BC::Expression>("velocity.bc", value.velocity_bc, AMREX_SPACEDIM);
    pp.select_default<BC::Constant,BC::Expression>("temperature.bc", value.temperature_bc, 1);
    pp.select_default<BC::Constant,BC::Expression>("mass_fraction.bc", value.mass_fraction_bc, value.nspecies);
    pp.select_default<BC::Constant,BC::Expression>("pressure.bc", value.pressure_bc, 1);
    pp.select_default<BC::Constant::ZeroNeumann, BC::Constant,BC::Expression>("eta.bc", value.eta_bc, 1);
    pp.select_default<BC::Constant::ZeroNeumann, BC::Constant,BC::Expression>("xi.bc", value.xi_bc, AMREX_SPACEDIM);

    pp.select_default<IC::Constant,IC::Expression>("velocity.ic", value.velocity_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("temperature.ic", value.temperature_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("mass_fraction.ic", value.mass_fraction_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("pressure.ic", value.pressure_ic, value.geom);
    if (pp.contains("eta.ic.type"))
        pp.select_default<IC::Constant,IC::Expression>("eta.ic", value.eta_ic, value.geom);
    else
        value.eta_ic = nullptr;
    pp.select_default<IC::Expression::X,IC::Constant,IC::Expression>("xi.ic", value.xi_ic, value.geom);

    value.AddField<Set::Scalar,Set::HC::Cell>(value.velocity_mf,        value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity",          true,  true, {"x","y"});
    value.AddField<Set::Scalar,Set::HC::Cell>(value.velocity_old_mf,    value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity_old",      false, true, {"x","y"});
    value.AddField<Set::Scalar,Set::HC::Cell>(value.temperature_mf,       value.temperature_bc,   1,              nghost, "temperature",       true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.temperature_old_mf,   value.temperature_bc,   1,              nghost, "temperature_old",   false, true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.mass_fraction_mf,     value.mass_fraction_bc, value.nspecies, nghost, "mass_fraction",     true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.mass_fraction_old_mf, value.mass_fraction_bc, value.nspecies, nghost, "mass_fraction_old", false, true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.eta_mf,               value.eta_bc,           1,              nghost, "eta",               true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.eta_old_mf,           value.eta_bc,           1,              nghost, "eta_old",           false, true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.xi_mf,                value.xi_bc,            AMREX_SPACEDIM, nghost, "xi",                true,  true, {"x","y"});
    value.AddField<Set::Scalar,Set::HC::Cell>(value.xi_old_mf,            value.xi_bc,            AMREX_SPACEDIM, nghost, "xi_old",            false, true, {"x","y"});

    value.AddField<Set::Scalar,Set::HC::Cell>(value.density_mf,             &value.bc_nothing, 1,              1,      "density",             true,  false);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.pressure_mf,            value.pressure_bc, 1,              nghost, "pressure",            true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.pressure_correction_mf, &value.bc_nothing, 1,              nghost, "pressure_correction", true, false);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.mole_fraction_mf, &value.bc_nothing, value.nspecies, 1, "mole_fraction", value.diagnostics_extended_fields, false);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.cell_weighted_solid_deviatoric_cauchy_stress_mf,         &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "cell_weighted_solid_deviatoric_cauchy_stress", true, false, {"_xx","_xy","_yx","_yy"});
    if (value.diagnostics_extended_fields)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(value.momentum_mf, &value.bc_nothing, AMREX_SPACEDIM, 1, "momentum", true, false, {"x","y"});
        value.AddField<Set::Scalar,Set::HC::Cell>(value.energy_mf, &value.bc_nothing, 1, 1, "energy", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.vorticity_mf, &value.bc_nothing, 1, 1, "vorticity", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.solid_weight_mf, &value.bc_nothing, 1, 1, "solid_weight", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.deformation_gradient_mf, &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "F", true, false, {"_xx","_xy","_yx","_yy"});
        value.AddField<Set::Scalar,Set::HC::Cell>(value.solid_first_piola_kirchhoff_stress_mf, &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "solid_first_piola_kirchhoff_stress", true, false, {"_xx","_xy","_yx","_yy"});
        value.AddField<Set::Scalar,Set::HC::Cell>(value.cell_total_cauchy_stress_mf, &value.bc_nothing, AMREX_SPACEDIM * AMREX_SPACEDIM, 1, "cell_total_cauchy_stress", true, false, {"_xx","_xy","_yx","_yy"});
    }

    bool allow_unused;
    pp.query_default("allow_unused", allow_unused, false);
    if (!allow_unused && pp.AnyUnusedInputs(true, false))
    {
        Util::Warning(INFO, "The following inputs were specified but not used:");
        pp.AllUnusedInputs();
        Util::Exception(INFO, "Aborting. Specify 'allow_unused=True` to ignore this error.");
    }
}

void
LowMach::SmoothReferenceMapForStress(int lev, const amrex::MultiFab& eta_mf, amrex::MultiFab& xi_mf)
{
    const int smooth_sweeps = reference_map_stress_smoothing_sweeps < 0 ? 0 : reference_map_stress_smoothing_sweeps;
    const Set::Scalar alpha = reference_map_stress_smoothing_alpha;
    if (smooth_sweeps == 0 || alpha == 0.0) return;

    const Set::Scalar stress_eta_min = Util::Clamp(reference_map_eta_extension, 0.0, 1.0);
    const Set::Scalar eta_core =
        Util::Max(Util::Clamp(reference_map_eta_core, 0.0, 1.0), 1.0e-12);
    const Set::Scalar smoothing_power = Util::Max(reference_map_smoothing_power, 0.0);
    const int xi_ngrow = xi_mf.nGrow();
    amrex::Geometry const geom_lev = geom[lev];
    amrex::Box const domain = geom[lev].Domain();
    const amrex::Dim3 lo = amrex::lbound(domain);
    const amrex::Dim3 hi = amrex::ubound(domain);

    for (int sweep = 0; sweep < smooth_sweeps; ++sweep)
    {
        xi_mf.FillBoundary(geom[lev].periodicity());

        amrex::MultiFab xi_old(xi_mf.boxArray(), xi_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        xi_old.FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(xi_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& eta = eta_mf.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_in = xi_old.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
                if (eta_val <= stress_eta_min) return;

                const Set::Scalar grade = 1.0 - Util::SmootherStep(eta_val / eta_core);
                if (grade == 0.0) return;
                const Set::Scalar local_alpha = alpha * std::pow(grade, smoothing_power);
                if (local_alpha == 0.0) return;

                Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar q = xi_in(i,j,k,n) - pos(n);
                    Set::Scalar lap = 0.0;

                    if (i > lo.x)
                    {
                        Set::Vector npos = Set::Position(i - 1, j, k, geom_lev, amrex::IndexType::TheCellType());
                        lap += xi_in(i-1,j,k,n) - npos(n) - q;
                    }
                    if (i < hi.x)
                    {
                        Set::Vector npos = Set::Position(i + 1, j, k, geom_lev, amrex::IndexType::TheCellType());
                        lap += xi_in(i+1,j,k,n) - npos(n) - q;
                    }
#if AMREX_SPACEDIM >= 2
                    if (j > lo.y)
                    {
                        Set::Vector npos = Set::Position(i, j - 1, k, geom_lev, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j-1,k,n) - npos(n) - q;
                    }
                    if (j < hi.y)
                    {
                        Set::Vector npos = Set::Position(i, j + 1, k, geom_lev, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j+1,k,n) - npos(n) - q;
                    }
#endif
#if AMREX_SPACEDIM == 3
                    if (k > lo.z)
                    {
                        Set::Vector npos = Set::Position(i, j, k - 1, geom_lev, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j,k-1,n) - npos(n) - q;
                    }
                    if (k < hi.z)
                    {
                        Set::Vector npos = Set::Position(i, j, k + 1, geom_lev, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j,k+1,n) - npos(n) - q;
                    }
#endif
                    xi(i,j,k,n) = pos(n) + q + local_alpha * lap;
                }
            });
        }
    }

    xi_mf.FillBoundary(geom[lev].periodicity());
}

void
LowMach::UpdateSolidStress(int lev,
                           const amrex::MultiFab& u_mf,
                           const amrex::MultiFab& T_mf,
                           const amrex::MultiFab& eta_mf,
                           const amrex::MultiFab& xi_mf,
                           bool write_diagnostics)
{
    BL_PROFILE("Integrator::LowMach::UpdateSolidStress");

    cell_weighted_solid_deviatoric_cauchy_stress_mf[lev]->setVal(0.0, 0, cell_weighted_solid_deviatoric_cauchy_stress_mf[lev]->nComp(), cell_weighted_solid_deviatoric_cauchy_stress_mf[lev]->nGrow());
    if (write_diagnostics)
    {
        deformation_gradient_mf[lev]->setVal(0.0);
        solid_first_piola_kirchhoff_stress_mf[lev]->setVal(0.0);
        cell_total_cauchy_stress_mf[lev]->setVal(0.0);
        solid_weight_mf[lev]->setVal(0.0);
    }

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    const Model::Solid::Finite::NeoHookean solid_model = finite_solid_model;
    const bool solid_enabled = finite_solid_enabled;
    const Set::Scalar eta_threshold = Util::Clamp(finite_solid_eta_threshold, 0.0, 1.0);
    const Set::Scalar stress_eta_min = Util::Clamp(reference_map_eta_extension, 0.0, 1.0);
    const Set::Scalar det_floor = Util::Max(small, finite_solid_J_floor);
    const Set::Scalar solid_viscosity = finite_solid_viscosity;
    const Set::Scalar solid_bulk_viscosity = finite_solid_bulk_viscosity;
    const Set::Scalar solid_interface_viscosity = finite_solid_interface_viscosity;
    const Set::Scalar p_scale = pressure_scale;
    const bool fluid_viscous = include_viscosity;

    amrex::MultiFab xi_stress_mf;
    const amrex::MultiFab* xi_for_stress_mf = &xi_mf;
    if (reference_map_stress_smoothing_sweeps > 0 && reference_map_stress_smoothing_alpha != 0.0)
    {
        xi_stress_mf.define(xi_mf.boxArray(), xi_mf.DistributionMap(), AMREX_SPACEDIM, xi_mf.nGrow());
        amrex::MultiFab::Copy(xi_stress_mf, xi_mf, 0, 0, AMREX_SPACEDIM, xi_mf.nGrow());
        SmoothReferenceMapForStress(lev, eta_mf, xi_stress_mf);
        xi_for_stress_mf = &xi_stress_mf;
    }

    for (amrex::MFIter mfi(*xi_for_stress_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        bx &= domain;
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.array(mfi);
        Set::Patch<const Set::Scalar> xi = xi_for_stress_mf->array(mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> cell_weighted_solid_deviatoric_cauchy_stress = cell_weighted_solid_deviatoric_cauchy_stress_mf.Patch(lev,mfi);
        amrex::Array4<Set::Scalar> F_field;
        amrex::Array4<Set::Scalar> solid_first_piola_kirchhoff_stress;
        amrex::Array4<Set::Scalar> cell_total_cauchy_stress;
        amrex::Array4<Set::Scalar> solid_weight_field;
        if (write_diagnostics)
        {
            F_field = deformation_gradient_mf[lev]->array(mfi);
            solid_first_piola_kirchhoff_stress = solid_first_piola_kirchhoff_stress_mf[lev]->array(mfi);
            cell_total_cauchy_stress = cell_total_cauchy_stress_mf[lev]->array(mfi);
            solid_weight_field = solid_weight_mf[lev]->array(mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix F = Set::Matrix::Identity();
            Set::Matrix P = Set::Matrix::Zero();
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Scalar div_u = grad_u.trace();

            Set::Matrix newtonian_viscous_stress = Set::Matrix::Zero();
            if (write_diagnostics && fluid_viscous)
            {
                Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
                newtonian_viscous_stress += mu * (grad_u + grad_u.transpose());
            }

            Set::Matrix solid_sigma = Set::Matrix::Zero();

            Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
            Set::Scalar solid_weight = 0.0;
            if (eta_val > eta_threshold)
            {
                Set::Scalar denom = Util::Max(1.0 - eta_threshold, 1.0e-12);
                solid_weight = Util::SmootherStep((eta_val - eta_threshold) / denom);
            }

            if (solid_enabled && eta_val > stress_eta_min)
            {
                Set::Matrix grad_xi = Numeric::Gradient(xi, i, j, k, DX, sten);
                Set::Scalar det_grad_xi = grad_xi.determinant();
                bool valid = Util::Abs(det_grad_xi) > det_floor;
                if (valid)
                {
#if AMREX_SPACEDIM == 2
                    Set::Scalar inv_det = 1.0 / det_grad_xi;
                    F(0,0) =  grad_xi(1,1) * inv_det;
                    F(0,1) = -grad_xi(0,1) * inv_det;
                    F(1,0) = -grad_xi(1,0) * inv_det;
                    F(1,1) =  grad_xi(0,0) * inv_det;
#else
                    F = grad_xi.inverse();
#endif
                }

                Set::Scalar J = F.determinant();
                if (valid && Util::Abs(J) > det_floor)
                {
                    P = solid_model.DW(F);
                    solid_sigma = (P * F.transpose()) / J;
                }

                if (solid_viscosity != 0.0 || solid_bulk_viscosity != 0.0)
                {
                    solid_sigma += solid_viscosity * (grad_u + grad_u.transpose()) +
                                   solid_bulk_viscosity * div_u * Set::Matrix::Identity();
                }
            }

            Set::Matrix solid_sigma_dev = solid_sigma -
                                          (solid_sigma.trace() / Set::Scalar(AMREX_SPACEDIM)) *
                                              Set::Matrix::Identity();
            Set::Matrix strain_rate_dev = grad_u + grad_u.transpose();
            strain_rate_dev -= (strain_rate_dev.trace() / Set::Scalar(AMREX_SPACEDIM)) *
                               Set::Matrix::Identity();
            Set::Matrix interface_damping_stress =
                solid_interface_viscosity * 4.0 * eta_val * (1.0 - eta_val) * strain_rate_dev;
            Set::Matrix weighted_solid_deviatoric_cauchy_stress =
                solid_weight * solid_sigma_dev + interface_damping_stress;

            cell_weighted_solid_deviatoric_cauchy_stress(i,j,k,0) = weighted_solid_deviatoric_cauchy_stress(0,0);
            cell_weighted_solid_deviatoric_cauchy_stress(i,j,k,1) = weighted_solid_deviatoric_cauchy_stress(0,1);
            cell_weighted_solid_deviatoric_cauchy_stress(i,j,k,2) = weighted_solid_deviatoric_cauchy_stress(1,0);
            cell_weighted_solid_deviatoric_cauchy_stress(i,j,k,3) = weighted_solid_deviatoric_cauchy_stress(1,1);
            if (write_diagnostics)
            {
                Set::Matrix total_cauchy_stress = -p_scale * pressure(i,j,k) * Set::Matrix::Identity() +
                                                  newtonian_viscous_stress +
                                                  solid_weight * solid_sigma +
                                                  interface_damping_stress;
                F_field(i,j,k,0) = F(0,0);
                F_field(i,j,k,1) = F(0,1);
                F_field(i,j,k,2) = F(1,0);
                F_field(i,j,k,3) = F(1,1);
                solid_first_piola_kirchhoff_stress(i,j,k,0) = P(0,0);
                solid_first_piola_kirchhoff_stress(i,j,k,1) = P(0,1);
                solid_first_piola_kirchhoff_stress(i,j,k,2) = P(1,0);
                solid_first_piola_kirchhoff_stress(i,j,k,3) = P(1,1);
                cell_total_cauchy_stress(i,j,k,0) = total_cauchy_stress(0,0);
                cell_total_cauchy_stress(i,j,k,1) = total_cauchy_stress(0,1);
                cell_total_cauchy_stress(i,j,k,2) = total_cauchy_stress(1,0);
                cell_total_cauchy_stress(i,j,k,3) = total_cauchy_stress(1,1);
                solid_weight_field(i,j,k) = solid_weight;
            }
        });
    }

    cell_weighted_solid_deviatoric_cauchy_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (write_diagnostics)
    {
        solid_weight_mf[lev]->FillBoundary(geom[lev].periodicity());
        deformation_gradient_mf[lev]->FillBoundary(geom[lev].periodicity());
        solid_first_piola_kirchhoff_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
        cell_total_cauchy_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
}

void
LowMach::UpdateThermodynamics(int lev, const amrex::MultiFab& T_mf, const amrex::MultiFab& Y_mf)
{
    amrex::Box domain = geom[lev].Domain();

    density_mf[lev]->setVal(0.0);
    mole_fraction_mf[lev]->setVal(0.0);

    for (amrex::MFIter mfi(Y_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        bx &= domain;
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> Y = Y_mf.array(mfi);
        Set::Patch<Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar p_floor = pressure_floor;
        const Set::Scalar p0 = thermodynamic_pressure;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar moles = 0.0;
            for (int n = 0; n < gas.nspecies; ++n) moles += Y(i,j,k,n) / gas.MW[n];
            if (!(moles > 0.0)) moles = 1.0 / gas.MW[0];
            for (int n = 0; n < gas.nspecies; ++n) X(i,j,k,n) = (Y(i,j,k,n) / gas.MW[n]) / moles;

            Set::Scalar p = Util::Max(p0, p_floor);
            Set::Scalar density = p / (gas.R(X, i, j, k) * T(i,j,k));
            density = Util::Max(density, rho_floor);
            rho(i,j,k) = density;
        });
    }

    density_mf[lev]->FillBoundary(geom[lev].periodicity());
    mole_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
}

void
LowMach::UpdateDerivedDiagnostics(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf)
{
    if (!diagnostics_extended_fields) return;

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    momentum_mf[lev]->setVal(0.0);
    energy_mf[lev]->setVal(0.0);
    vorticity_mf[lev]->setVal(0.0);

    for (amrex::MFIter mfi(u_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        bx &= domain;
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> M = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> E = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> omega = vorticity_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Set::Scalar density = rho(i,j,k);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                M(i,j,k,d) = density * u(i,j,k,d);
            E(i,j,k) = gas.ComputeE(density, M(i,j,k,0), M(i,j,k,1), T(i,j,k), X, i, j, k);

            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            omega(i,j,k) = grad_u(1,0) - grad_u(0,1);
        });
    }

    momentum_mf[lev]->FillBoundary(geom[lev].periodicity());
    energy_mf[lev]->FillBoundary(geom[lev].periodicity());
    vorticity_mf[lev]->FillBoundary(geom[lev].periodicity());
}


void
LowMach::ProjectVelocity(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ProjectVelocity");
    if (!projection_enabled || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    amrex::Vector<amrex::Geometry> proj_geom(nlev);
    amrex::Vector<amrex::BoxArray> proj_grids(nlev);
    amrex::Vector<amrex::DistributionMapping> proj_dmap(nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> phi(nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhs(nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> beta_cc(nlev);
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> beta_face(nlev);
    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>> beta_face_const(nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        proj_geom[lev] = geom[lev];
        proj_grids[lev] = velocity_mf[lev]->boxArray();
        proj_dmap[lev] = velocity_mf[lev]->DistributionMap();

        velocity_bc->define(geom[lev]);
        velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, time, 0);
        velocity_mf[lev]->FillBoundary(geom[lev].periodicity());

        UpdateThermodynamics(lev, *temperature_mf[lev], *mass_fraction_mf[lev]);

        phi[lev].reset(new amrex::MultiFab(proj_grids[lev], proj_dmap[lev], 1, pressure_correction_mf[lev]->nGrow()));
        rhs[lev].reset(new amrex::MultiFab(proj_grids[lev], proj_dmap[lev], 1, 0));
        beta_cc[lev].reset(new amrex::MultiFab(proj_grids[lev], proj_dmap[lev], 1, 1));
        phi[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);
        beta_cc[lev]->setVal(0.0);
    }

    for (int lev = 0; lev < nlev; ++lev)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();
        amrex::MultiFab& u_mf = *velocity_mf[lev];
        amrex::MultiFab& rho_mf = *density_mf[lev];
        const Set::Scalar inv_dt = 1.0 / dt;
        const Set::Scalar rho_floor = density_floor;

        for (amrex::MFIter mfi(*rhs[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
            Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
            Set::Patch<Set::Scalar> rhs_arr = rhs[lev]->array(mfi);
            Set::Patch<Set::Scalar> beta = beta_cc[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, domain);
                Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
                Set::Scalar div_u = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) div_u += grad_u(d,d);
                rhs_arr(i,j,k) = div_u * inv_dt;
                beta(i,j,k) = 1.0 / Util::Max(rho(i,j,k), rho_floor);
            });
        }

        beta_cc[lev]->FillBoundary(geom[lev].periodicity());
    }

    BC::Constant::ZeroNeumann beta_bc(1);
    for (int lev = 0; lev < nlev; ++lev)
    {
        beta_bc.define(geom[lev]);
        if (lev == 0)
        {
            beta_bc.FillBoundary(*beta_cc[lev], 0, 1, time, 0);
            beta_cc[lev]->FillBoundary(geom[lev].periodicity());
        }
        else
        {
            amrex::Vector<amrex::MultiFab*> cmf{beta_cc[lev - 1].get()};
            amrex::Vector<amrex::MultiFab*> fmf{beta_cc[lev].get()};
            amrex::Vector<amrex::Real> ctime{time};
            amrex::Vector<amrex::Real> ftime{time};
            amrex::Vector<amrex::BCRec> bcs(1, beta_bc.GetBCRec());
            amrex::FillPatchTwoLevels(*beta_cc[lev], time, cmf, ctime, fmf, ftime,
                                       0, 0, 1, geom[lev - 1], geom[lev],
                                       beta_bc, 0, beta_bc, 0,
                                       refRatio(lev - 1), &amrex::cell_cons_interp, bcs, 0);
        }
    }

    for (int lev = 0; lev < nlev; ++lev)
    {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> beta_face_ptr;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_ba = proj_grids[lev];
            face_ba.surroundingNodes(d);
            beta_face[lev][d].define(face_ba, proj_dmap[lev], 1, 0);
            beta_face_ptr[d] = &beta_face[lev][d];
            beta_face_const[lev][d] = &beta_face[lev][d];
        }
        amrex::average_cellcenter_to_face(beta_face_ptr, *beta_cc[lev], geom[lev], 1, true, 0);
    }

    amrex::LPInfo info;
    amrex::MLABecLaplacian mlabec(proj_geom, proj_grids, proj_dmap, info);
    mlabec.setMaxOrder(2);
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> lobc;
    amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> hibc;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        lobc[d] = geom[0].isPeriodic(d) ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
        hibc[d] = geom[0].isPeriodic(d) ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Neumann;
    }
    mlabec.setDomainBC(lobc, hibc);
    mlabec.setScalars(0.0, -1.0);
    for (int lev = 0; lev < nlev; ++lev)
    {
        mlabec.setLevelBC(lev, nullptr);
        mlabec.setACoeffs(lev, 0.0);
        mlabec.setBCoeffs(lev, beta_face_const[lev]);
    }

    amrex::Vector<amrex::MultiFab*> phi_ptr(nlev);
    amrex::Vector<amrex::MultiFab const*> rhs_ptr(nlev);
    for (int lev = 0; lev < nlev; ++lev)
    {
        phi_ptr[lev] = phi[lev].get();
        rhs_ptr[lev] = rhs[lev].get();
    }

    amrex::MLMG mlmg(mlabec);
    mlmg.setVerbose(projection_verbose);
    mlmg.solve(phi_ptr, rhs_ptr, projection_tol_rel, projection_tol_abs);

    BC::Constant::ZeroNeumann phi_bc(1);
    for (int lev = 0; lev < nlev; ++lev)
    {
        phi_bc.define(geom[lev]);
        if (lev == 0)
        {
            phi_bc.FillBoundary(*phi[lev], 0, 1, time, 0);
            phi[lev]->FillBoundary(geom[lev].periodicity());
        }
        else
        {
            amrex::Vector<amrex::MultiFab*> cmf{phi[lev - 1].get()};
            amrex::Vector<amrex::MultiFab*> fmf{phi[lev].get()};
            amrex::Vector<amrex::Real> ctime{time};
            amrex::Vector<amrex::Real> ftime{time};
            amrex::Vector<amrex::BCRec> bcs(1, phi_bc.GetBCRec());
            amrex::FillPatchTwoLevels(*phi[lev], time, cmf, ctime, fmf, ftime,
                                       0, 0, 1, geom[lev - 1], geom[lev],
                                       phi_bc, 0, phi_bc, 0,
                                       refRatio(lev - 1), &amrex::cell_cons_interp, bcs, 0);
        }
    }

    for (int lev = 0; lev < nlev; ++lev)
    {
        amrex::MultiFab::Copy(*pressure_correction_mf[lev], *phi[lev], 0, 0, 1, pressure_correction_mf[lev]->nGrow());
        pressure_correction_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
    for (int lev = 0; lev < nlev; ++lev)
    {
        amrex::MultiFab& u_mf = *velocity_mf[lev];
        amrex::MultiFab& p_mf = *pressure_mf[lev];
        amrex::MultiFab& rho_mf = *density_mf[lev];
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar p_floor = pressure_floor;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar p_scale_inv = 1.0 / p_scale;
        const bool update_pressure = projection_update_pressure;

        for (amrex::MFIter mfi(u_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> u = u_mf.array(mfi);
            Set::Patch<Set::Scalar> p = p_mf.array(mfi);
            Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
            Set::Patch<const Set::Scalar> phi_arr = pressure_correction_mf[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, domain);
                Set::Vector grad_phi = Numeric::Gradient(phi_arr, i, j, k, 0, DX, sten);
                const Set::Scalar beta = 1.0 / Util::Max(rho(i,j,k), rho_floor);
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    u(i,j,k,d) -= dt * beta * grad_phi(d);
                if (update_pressure) p(i,j,k) = Util::Max(p(i,j,k) + p_scale_inv * phi_arr(i,j,k), p_floor);
            });
        }

        velocity_bc->define(geom[lev]);
        velocity_bc->FillBoundary(u_mf, 0, AMREX_SPACEDIM, time, 0);
        u_mf.FillBoundary(geom[lev].periodicity());
        pressure_bc->define(geom[lev]);
        pressure_bc->FillBoundary(p_mf, 0, 1, time, 0);
        p_mf.FillBoundary(geom[lev].periodicity());
    }

    for (int lev = nlev - 2; lev >= 0; --lev)
        amrex::average_down(*velocity_mf[lev + 1], *velocity_mf[lev],
                            geom[lev + 1], geom[lev], 0, AMREX_SPACEDIM, refRatio(lev));
}

void
LowMach::RebuildReferenceMapOutsideEta(int lev, const amrex::MultiFab& eta_stage_mf, amrex::MultiFab& xi_stage_mf, Set::Scalar time)
{
    const Set::Scalar core_eta =
        Util::Max(Util::Clamp(reference_map_eta_core, 0.0, 1.0), 1.0e-12);
    const Set::Scalar reconstruction_alpha =
        Util::Clamp(reference_map_reconstruction_alpha, 0.0, 1.0);
    const Set::Scalar reconstruction_power =
        Util::Max(reference_map_reconstruction_power, 0.0);
    const Set::Scalar affine_tolerance =
        Util::Max(reference_map_affine_tolerance, 1.0e-12);
    const Set::Scalar smoothing_alpha =
        Util::Clamp(reference_map_smoothing_alpha, 0.0, 1.0);
    const Set::Scalar smoothing_power = Util::Max(reference_map_smoothing_power, 0.0);
    const int extrap_sweeps = reference_map_extrapolation_sweeps < 0 ? 0 : reference_map_extrapolation_sweeps;
    const int smooth_sweeps = reference_map_smoothing_sweeps < 0 ? 0 : reference_map_smoothing_sweeps;
    const int xi_ngrow = xi_stage_mf.nGrow();
    const amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM> DX = geom[lev].CellSizeArray();

    amrex::MultiFab eta_work(eta_stage_mf.boxArray(), eta_stage_mf.DistributionMap(), 1, eta_stage_mf.nGrow());
    amrex::MultiFab::Copy(eta_work, eta_stage_mf, 0, 0, 1, eta_stage_mf.nGrow());
    eta_work.FillBoundary(geom[lev].periodicity());

    auto fill_xi_boundary = [&]()
    {
        xi_bc->define(geom[lev]);
        xi_bc->FillBoundary(xi_stage_mf, 0, AMREX_SPACEDIM, time, 0);
        xi_stage_mf.FillBoundary(geom[lev].periodicity());
    };

    fill_xi_boundary();

    for (int sweep = 0; sweep < extrap_sweeps; ++sweep)
    {
        amrex::MultiFab xi_old(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_stage_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        xi_old.FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& eta = eta_work.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_in = xi_old.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
                const Set::Scalar grade =
                    1.0 - Util::SmootherStep(eta_val / core_eta);
                if (grade == 0.0) return;
                const Set::Scalar repair =
                    reconstruction_alpha * std::pow(grade, reconstruction_power);
                if (repair == 0.0) return;

                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar sum = 0.0;
                    Set::Scalar weight = 0.0;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        for (int sign = -1; sign <= 1; sign += 2)
                        {
                            const int ii = i + (d == 0 ? sign : 0);
                            const int jj = j + (d == 1 ? sign : 0);
                            const int kk = k + (d == 2 ? sign : 0);
                            const Set::Scalar eta_near = Util::Clamp(eta(ii,jj,kk), 0.0, 1.0);
                            if (eta_near <= eta_val) continue;

                            const int iii = i + (d == 0 ? 2 * sign : 0);
                            const int jjj = j + (d == 1 ? 2 * sign : 0);
                            const int kkk = k + (d == 2 ? 2 * sign : 0);
                            const Set::Scalar eta_far = Util::Clamp(eta(iii,jjj,kkk), 0.0, 1.0);
                            const Set::Scalar nearest =
                                xi_in(ii,jj,kk,n) - (n == d ? sign * DX[d] : 0.0);
                            Set::Scalar candidate = nearest;
                            if (eta_far > eta_near)
                            {
                                const Set::Scalar inward_slope =
                                    xi_in(ii,jj,kk,n) - xi_in(iii,jjj,kkk,n);
                                const Set::Scalar linear =
                                    xi_in(ii,jj,kk,n) + inward_slope;
                                const Set::Scalar eta_weight =
                                    Util::Clamp(eta_near / core_eta, 0.0, 1.0);
                                Set::Scalar slope_weight = eta_weight;

                                const int iiii = i + (d == 0 ? 3 * sign : 0);
                                const int jjjj = j + (d == 1 ? 3 * sign : 0);
                                const int kkkk = k + (d == 2 ? 3 * sign : 0);
                                const Set::Scalar eta_farther =
                                    Util::Clamp(eta(iiii,jjjj,kkkk), 0.0, 1.0);
                                if (eta_farther >= eta_far)
                                {
                                    Set::Scalar curvature_sq = 0.0;
                                    Set::Scalar slope_sq = 0.0;
                                    for (int m = 0; m < AMREX_SPACEDIM; ++m)
                                    {
                                        const Set::Scalar inward_vector_slope =
                                            xi_in(ii,jj,kk,m) - xi_in(iii,jjj,kkk,m);
                                        const Set::Scalar farther_vector_slope =
                                            xi_in(iii,jjj,kkk,m) - xi_in(iiii,jjjj,kkkk,m);
                                        const Set::Scalar curvature =
                                            inward_vector_slope - farther_vector_slope;
                                        curvature_sq += curvature * curvature;
                                        slope_sq += inward_vector_slope * inward_vector_slope +
                                                    farther_vector_slope * farther_vector_slope;
                                    }
                                    const Set::Scalar relative_curvature =
                                        std::sqrt(curvature_sq) /
                                        (std::sqrt(slope_sq) + 1.0e-12);
                                    const Set::Scalar affine_weight = 1.0 - Util::SmootherStep(
                                        relative_curvature / affine_tolerance);
                                    slope_weight = eta_weight +
                                        (1.0 - eta_weight) * affine_weight;
                                }
                                candidate = nearest + slope_weight * (linear - nearest);
                            }
                            const Set::Scalar candidate_weight = eta_near - eta_val;
                            sum += candidate_weight * candidate;
                            weight += candidate_weight;
                        }
                    }

                    if (weight > 1.0e-14)
                    {
                        const Set::Scalar reconstructed = sum / weight;
                        xi(i,j,k,n) = (1.0 - repair) * xi_in(i,j,k,n) + repair * reconstructed;
                    }
                }
            });
        }
        amrex::Gpu::streamSynchronize();
        fill_xi_boundary();
    }

    for (int sweep = 0; sweep < smooth_sweeps; ++sweep)
    {
        amrex::MultiFab xi_old(xi_stage_mf.boxArray(), xi_stage_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_stage_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        xi_old.FillBoundary(geom[lev].periodicity());

        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& eta = eta_work.const_array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_in = xi_old.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
                const Set::Scalar grade = 1.0 - Util::SmootherStep(eta_val / core_eta);
                if (grade == 0.0) return;
                const Set::Scalar relax = smoothing_alpha * std::pow(grade, smoothing_power);
                if (relax == 0.0) return;

                Set::Scalar eta_support = eta_val +
                    Util::Clamp(eta(i-1,j,k), 0.0, 1.0) +
                    Util::Clamp(eta(i+1,j,k), 0.0, 1.0);
#if AMREX_SPACEDIM >= 2
                eta_support += Util::Clamp(eta(i,j-1,k), 0.0, 1.0) +
                               Util::Clamp(eta(i,j+1,k), 0.0, 1.0);
#endif
#if AMREX_SPACEDIM == 3
                eta_support += Util::Clamp(eta(i,j,k-1), 0.0, 1.0) +
                               Util::Clamp(eta(i,j,k+1), 0.0, 1.0);
#endif
                if (eta_support == 0.0) return;

                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar sum = xi_in(i-1,j,k,n) + xi_in(i+1,j,k,n);
                    Set::Scalar count = 2.0;
#if AMREX_SPACEDIM >= 2
                    sum += xi_in(i,j-1,k,n) + xi_in(i,j+1,k,n);
                    count += 2.0;
#endif
#if AMREX_SPACEDIM == 3
                    sum += xi_in(i,j,k-1,n) + xi_in(i,j,k+1,n);
                    count += 2.0;
#endif
                    xi(i,j,k,n) = (1.0 - relax) * xi_in(i,j,k,n) + relax * sum / count;
                }
            });
        }
        amrex::Gpu::streamSynchronize();
        fill_xi_boundary();
    }

    fill_xi_boundary();
}

void
LowMach::Initialize(int lev)
{
    BL_PROFILE("Integrator::LowMach::Initialize");

    velocity_mf[lev]->setVal(0.0, 0, velocity_mf[lev]->nComp(), velocity_mf[lev]->nGrow());
    velocity_old_mf[lev]->setVal(0.0, 0, velocity_old_mf[lev]->nComp(), velocity_old_mf[lev]->nGrow());
    temperature_mf[lev]->setVal(0.0, 0, temperature_mf[lev]->nComp(), temperature_mf[lev]->nGrow());
    temperature_old_mf[lev]->setVal(0.0, 0, temperature_old_mf[lev]->nComp(), temperature_old_mf[lev]->nGrow());
    mass_fraction_mf[lev]->setVal(0.0, 0, mass_fraction_mf[lev]->nComp(), mass_fraction_mf[lev]->nGrow());
    mass_fraction_old_mf[lev]->setVal(0.0, 0, mass_fraction_old_mf[lev]->nComp(), mass_fraction_old_mf[lev]->nGrow());
    eta_mf[lev]->setVal(0.0, 0, eta_mf[lev]->nComp(), eta_mf[lev]->nGrow());
    eta_old_mf[lev]->setVal(0.0, 0, eta_old_mf[lev]->nComp(), eta_old_mf[lev]->nGrow());
    xi_mf[lev]->setVal(0.0, 0, xi_mf[lev]->nComp(), xi_mf[lev]->nGrow());
    xi_old_mf[lev]->setVal(0.0, 0, xi_old_mf[lev]->nComp(), xi_old_mf[lev]->nGrow());
    pressure_mf[lev]->setVal(0.0, 0, pressure_mf[lev]->nComp(), pressure_mf[lev]->nGrow());

    velocity_ic->Initialize(lev, velocity_mf, 0.0);
    velocity_ic->Initialize(lev, velocity_old_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_old_mf, 0.0);
    mass_fraction_ic->Initialize(lev, mass_fraction_mf, 0.0);
    mass_fraction_ic->Initialize(lev, mass_fraction_old_mf, 0.0);
    if (eta_ic)
    {
        eta_ic->Initialize(lev, eta_mf, 0.0);
        eta_ic->Initialize(lev, eta_old_mf, 0.0);
    }
    else
    {
        eta_mf[lev]->setVal(eta_initial_value);
        eta_old_mf[lev]->setVal(eta_initial_value);
    }
    xi_ic->Initialize(lev, xi_mf, 0.0);
    xi_ic->Initialize(lev, xi_old_mf, 0.0);
    pressure_ic->Initialize(lev, pressure_mf, 0.0);
    pressure_correction_mf[lev]->setVal(0.0);

    velocity_bc->define(geom[lev]);
    temperature_bc->define(geom[lev]);
    mass_fraction_bc->define(geom[lev]);
    eta_bc->define(geom[lev]);
    xi_bc->define(geom[lev]);
    pressure_bc->define(geom[lev]);

    velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
    temperature_bc->FillBoundary(*temperature_mf[lev], 0, 1, 0.0, 0);
    temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
    mass_fraction_bc->FillBoundary(*mass_fraction_mf[lev], 0, nspecies, 0.0, 0);
    mass_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
    eta_bc->FillBoundary(*eta_mf[lev], 0, 1, 0.0, 0);
    eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    xi_bc->FillBoundary(*xi_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    xi_mf[lev]->FillBoundary(geom[lev].periodicity());

    velocity_bc->FillBoundary(*velocity_old_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    velocity_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    temperature_bc->FillBoundary(*temperature_old_mf[lev], 0, 1, 0.0, 0);
    temperature_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    mass_fraction_bc->FillBoundary(*mass_fraction_old_mf[lev], 0, nspecies, 0.0, 0);
    mass_fraction_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    eta_bc->FillBoundary(*eta_old_mf[lev], 0, 1, 0.0, 0);
    eta_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    xi_bc->FillBoundary(*xi_old_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    xi_old_mf[lev]->FillBoundary(geom[lev].periodicity());

    pressure_bc->FillBoundary(*pressure_mf[lev], 0, 1, 0.0, 0);
    pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
    RebuildReferenceMapOutsideEta(lev, *eta_mf[lev], *xi_mf[lev], 0.0);
    RebuildReferenceMapOutsideEta(lev, *eta_old_mf[lev], *xi_old_mf[lev], 0.0);
    UpdateThermodynamics(lev, *temperature_mf[lev], *mass_fraction_mf[lev]);
    UpdateSolidStress(lev, *velocity_mf[lev], *temperature_mf[lev], *eta_mf[lev], *xi_mf[lev], diagnostics_extended_fields);
    UpdateDerivedDiagnostics(lev, *velocity_mf[lev], *temperature_mf[lev]);
}

void
LowMach::RHS(int lev, Set::Scalar /*time*/,
             amrex::MultiFab& u_rhs_mf,
             amrex::MultiFab& T_rhs_mf,
             amrex::MultiFab& Y_rhs_mf,
             amrex::MultiFab& eta_rhs_mf,
             amrex::MultiFab& xi_rhs_mf,
             const amrex::MultiFab& u_mf,
             const amrex::MultiFab& T_mf,
             const amrex::MultiFab& Y_mf,
             const amrex::MultiFab& eta_mf,
             const amrex::MultiFab& xi_mf)
{
    UpdateThermodynamics(lev, T_mf, Y_mf);
    if (finite_solid_deviatoric_stress_divergence_sign != 0.0)
        UpdateSolidStress(lev, u_mf, T_mf, eta_mf, xi_mf);

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    const auto advect_op = advect;
    const Numeric::Advect::Options advective_options{Numeric::Advect::Form::Advective};

    u_rhs_mf.setVal(0.0, 0, AMREX_SPACEDIM, u_rhs_mf.nGrow());

    for (amrex::MFIter mfi(u_rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> cell_weighted_solid_deviatoric_cauchy_stress =
            cell_weighted_solid_deviatoric_cauchy_stress_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> u_rhs = u_rhs_mf.array(mfi);
        const bool viscous = include_viscosity;
        const Set::Vector gravity = g;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar rho_floor = density_floor;
        const bool pressure_predictor = projection_predictor_pressure_gradient;
        const Set::Scalar deviatoric_stress_divergence_sign = finite_solid_deviatoric_stress_divergence_sign;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
            Set::Scalar mu = viscous ? gas.dynamic_viscosity(T(i,j,k), X, i, j, k) : 0.0;

            Set::Vector grad_p_vec = Set::Vector::Zero();
            if (pressure_predictor)
                grad_p_vec = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
            Set::Vector div_sigma = Set::Vector::Zero();
            if (deviatoric_stress_divergence_sign != 0.0)
                div_sigma = Numeric::MatrixDivergence(cell_weighted_solid_deviatoric_cauchy_stress, i, j, k, 0, DX, sten);

            Set::Vector adv_vec = advect_op.Vector(u, u, i, j, k, 0, DX, advective_options, sten);
            Set::Vector lap_vec = Set::Vector::Zero();
            if (viscous)
                lap_vec = Numeric::VectorLaplacian(u, i, j, k, 0, DX);

            Set::Vector rhs_vec = adv_vec
                                  - (p_scale / density) * grad_p_vec
                                  + gravity
                                  + (mu / density) * lap_vec
                                  + (deviatoric_stress_divergence_sign / density) * div_sigma;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                u_rhs(i,j,k,d) = rhs_vec(d);
        });
    }

    for (amrex::MFIter mfi(T_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> Y = Y_mf.array(mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.array(mfi);
        Set::Patch<const Set::Scalar> xi = xi_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> Y_rhs = Y_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> eta_rhs = eta_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> xi_rhs = xi_rhs_mf.array(mfi);
        const int nsp = nspecies;
        const bool conductive = include_conduction;
        const bool advect_T = advect_temperature;
        const Set::Scalar rho_floor = density_floor;
        const bool eta_phase_enabled = eta_phase_field_enabled;
        const Set::Scalar eta_phase_epsilon = eta_phase_field_epsilon;
        const Set::Scalar eta_phase_mobility = eta_phase_field_mobility;
        const Set::Scalar eta_phase_counter_curvature = eta_phase_field_counter_curvature;
        const Set::Scalar eta_phase_band = eta_phase_field_band;
        const bool eta_phase_mapped_distance = eta_phase_field_mapped_distance;
        const Set::Scalar eta_phase_normal_coherence_threshold =
            eta_phase_field_normal_coherence_threshold;
        const Set::Scalar eta_phase_normal_coherence_power =
            eta_phase_field_normal_coherence_power;
        const Set::Scalar eta_phase_exterior_reinitialization =
            eta_phase_field_exterior_reinitialization;
        const Set::Scalar eta_phase_incoherent_diffusion =
            eta_phase_field_incoherent_diffusion;
        const Set::Scalar eta_phase_normal_regularization =
            eta_phase_field_normal_regularization /
            (eta_phase_mapped_distance ? 1.0 : Util::Max(eta_phase_field_epsilon, small));

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Vector vel = Set::Vector::Zero();
            vel(0) = u(i,j,k,0);
            vel(1) = u(i,j,k,1);
#if AMREX_SPACEDIM == 3
            vel(2) = u(i,j,k,2);
#endif
            Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);

            T_rhs(i,j,k) = 0.0;
            if (advect_T)
            {
                Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
                T_rhs(i,j,k) -= vel.dot(grad_T);
            }
            if (conductive)
            {
                Set::Scalar kappa = gas.thermal_conductivity(T(i,j,k), X, i, j, k);
                Set::Scalar cp = gas.cp_mass(T(i,j,k), X, i, j, k);
                if (kappa == kappa && cp == cp && cp > 0.0)
                    T_rhs(i,j,k) += kappa / (density * cp) * Numeric::Laplacian(T, i, j, k, 0, DX);
            }

            for (int n = 0; n < nsp; ++n)
            {
                Set::Vector grad_Y = Numeric::Gradient(Y, i, j, k, n, DX, sten);
                Y_rhs(i,j,k,n) = -vel.dot(grad_Y);
            }

            eta_rhs(i,j,k) = advect_op(eta, u, i, j, k, 0, DX, advective_options, sten);
            if (eta_phase_enabled)
                eta_rhs(i,j,k) += EtaPhaseFieldSource(eta, i, j, k, DX,
                                                                  eta_phase_epsilon,
                                                                  eta_phase_mobility,
                                                                  eta_phase_counter_curvature,
                                                                  eta_phase_band,
                                                                  eta_phase_normal_regularization,
                                                                  eta_phase_mapped_distance,
                                                                  eta_phase_normal_coherence_threshold,
                                                                  eta_phase_normal_coherence_power,
                                                                  eta_phase_exterior_reinitialization,
                                                                  eta_phase_incoherent_diffusion);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                xi_rhs(i,j,k,d) = advect_op(xi, u, i, j, k, d, DX, advective_options, sten);
        });
    }
    u_rhs_mf.FillBoundary(geom[lev].periodicity());
}

void
LowMach::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(velocity_old_mf[lev], velocity_mf[lev]);
    std::swap(temperature_old_mf[lev], temperature_mf[lev]);
    std::swap(mass_fraction_old_mf[lev], mass_fraction_mf[lev]);
    std::swap(eta_old_mf[lev], eta_mf[lev]);
    std::swap(xi_old_mf[lev], xi_mf[lev]);

    amrex::Vector<amrex::MultiFab> solution_new;
    solution_new.emplace_back(*velocity_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_new.emplace_back(*temperature_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_new.emplace_back(*mass_fraction_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);
    solution_new.emplace_back(*eta_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_new.emplace_back(*xi_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);

    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*velocity_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_old.emplace_back(*temperature_old_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_old.emplace_back(*mass_fraction_old_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);
    solution_old.emplace_back(*eta_old_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_old.emplace_back(*xi_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);

    velocity_bc->define(geom[lev]);
    temperature_bc->define(geom[lev]);
    mass_fraction_bc->define(geom[lev]);
    eta_bc->define(geom[lev]);
    xi_bc->define(geom[lev]);

    amrex::TimeIntegrator timeintegrator(solution_new, time);
    timeintegrator.set_rhs([&](amrex::Vector<amrex::MultiFab>& rhs_mf,
                               amrex::Vector<amrex::MultiFab>& state_mf,
                               const Set::Scalar rhs_time)
    {
        velocity_bc->FillBoundary(state_mf[0], 0, AMREX_SPACEDIM, rhs_time, 0);
        state_mf[0].FillBoundary(geom[lev].periodicity());
        temperature_bc->FillBoundary(state_mf[1], 0, 1, rhs_time, 0);
        state_mf[1].FillBoundary(geom[lev].periodicity());
        mass_fraction_bc->FillBoundary(state_mf[2], 0, nspecies, rhs_time, 0);
        state_mf[2].FillBoundary(geom[lev].periodicity());
        eta_bc->FillBoundary(state_mf[3], 0, 1, rhs_time, 0);
        state_mf[3].FillBoundary(geom[lev].periodicity());
        xi_bc->FillBoundary(state_mf[4], 0, AMREX_SPACEDIM, rhs_time, 0);
        state_mf[4].FillBoundary(geom[lev].periodicity());

        RebuildReferenceMapOutsideEta(lev, state_mf[3], state_mf[4], rhs_time);
        RHS(lev, rhs_time, rhs_mf[0], rhs_mf[1], rhs_mf[2], rhs_mf[3], rhs_mf[4],
            state_mf[0], state_mf[1], state_mf[2], state_mf[3], state_mf[4]);
    });

    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab>& stage_mf, Set::Scalar stage_time)
    {
        velocity_bc->FillBoundary(stage_mf[0], 0, AMREX_SPACEDIM, stage_time, 0);
        stage_mf[0].FillBoundary(geom[lev].periodicity());
        temperature_bc->FillBoundary(stage_mf[1], 0, 1, stage_time, 0);
        stage_mf[1].FillBoundary(geom[lev].periodicity());
        mass_fraction_bc->FillBoundary(stage_mf[2], 0, nspecies, stage_time, 0);
        stage_mf[2].FillBoundary(geom[lev].periodicity());
        eta_bc->FillBoundary(stage_mf[3], 0, 1, stage_time, 0);
        stage_mf[3].FillBoundary(geom[lev].periodicity());
        xi_bc->FillBoundary(stage_mf[4], 0, AMREX_SPACEDIM, stage_time, 0);
        stage_mf[4].FillBoundary(geom[lev].periodicity());

        RebuildReferenceMapOutsideEta(lev, stage_mf[3], stage_mf[4], stage_time);
    });

    timeintegrator.advance(solution_old, solution_new, time, dt);
    Set::Scalar new_time = time + dt;
    velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, new_time, 0);
    velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
    temperature_bc->FillBoundary(*temperature_mf[lev], 0, 1, new_time, 0);
    temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
    mass_fraction_bc->FillBoundary(*mass_fraction_mf[lev], 0, nspecies, new_time, 0);
    mass_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
    eta_bc->FillBoundary(*eta_mf[lev], 0, 1, new_time, 0);
    eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    xi_bc->FillBoundary(*xi_mf[lev], 0, AMREX_SPACEDIM, new_time, 0);
    xi_mf[lev]->FillBoundary(geom[lev].periodicity());

    RebuildReferenceMapOutsideEta(lev, *eta_mf[lev], *xi_mf[lev], new_time);
    eta_bc->FillBoundary(*eta_mf[lev], 0, 1, new_time, 0);
    eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    xi_bc->FillBoundary(*xi_mf[lev], 0, AMREX_SPACEDIM, new_time, 0);
    xi_mf[lev]->FillBoundary(geom[lev].periodicity());
    RebuildReferenceMapOutsideEta(lev, *eta_mf[lev], *xi_mf[lev], new_time);
}

void
LowMach::TimeStepBegin(Set::Scalar /*time*/, int /*iter*/)
{
    if (!dynamictimestep.on) return;

    Set::Scalar advmax = 0.0;
    Set::Scalar viscmax = 0.0;
    Set::Scalar elasticmax = 0.0;
    Set::Scalar phasefieldmax = 0.0;
    const bool explicit_solid_deviatoric_stress = finite_solid_enabled &&
                                       finite_solid_deviatoric_stress_divergence_sign != 0.0;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        UpdateThermodynamics(lev, *temperature_mf[lev], *mass_fraction_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dxmin = std::min(DX[0], DX[1]);
        if (eta_phase_field_enabled && eta_phase_field_mobility > 0.0 && eta_phase_field_epsilon > 0.0)
            phasefieldmax = std::max(phasefieldmax,
                                     eta_phase_field_mobility / dxmin
                                     + 0.5 * eta_phase_field_mobility * eta_phase_field_epsilon / (dxmin * dxmin));

        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                Set::Scalar speed = std::sqrt(u(i,j,k,0)*u(i,j,k,0) + u(i,j,k,1)*u(i,j,k,1));
                return {speed / dxmin};
            });
            ReduceTuple hv = reduce_data.value();
            advmax = std::max(advmax, amrex::get<0>(hv));
        }

        for (amrex::MFIter mfi(*temperature_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);
            const Set::Scalar rho_floor = density_floor;
            const bool viscous = include_viscosity;
            const bool elastic = explicit_solid_deviatoric_stress;
            const Set::Scalar eta_threshold = Util::Clamp(finite_solid_eta_threshold, 0.0, 1.0);
            const Set::Scalar mu_solid = finite_solid_model.mu;
            const Set::Scalar kappa_solid = finite_solid_model.kappa;
            const Set::Scalar solid_viscosity = finite_solid_viscosity;
            const Set::Scalar solid_bulk_viscosity = finite_solid_bulk_viscosity;
            const Set::Scalar solid_interface_viscosity = finite_solid_interface_viscosity;

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                Set::Scalar nu = 0.0;
                if (viscous)
                {
                    Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
                    nu = mu / Util::Max(rho(i,j,k), rho_floor);
                }
                Set::Scalar elastic_rate = 0.0;
                if (elastic && Util::Clamp(eta(i,j,k), 0.0, 1.0) > eta_threshold)
                {
                    Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
                    Set::Scalar wave_speed = std::sqrt(Util::Max(kappa_solid + (4.0 / 3.0) * mu_solid, mu_solid) / density);
                    elastic_rate = wave_speed / dxmin;
                    nu = Util::Max(nu, (solid_viscosity + solid_bulk_viscosity + solid_interface_viscosity) / density);
                }
                return {nu / (dxmin * dxmin), elastic_rate};
            });
            ReduceTuple hv = reduce_data.value();
            viscmax = std::max(viscmax, amrex::get<0>(hv));
            elasticmax = std::max(elasticmax, amrex::get<1>(hv));
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(advmax);
    amrex::ParallelDescriptor::ReduceRealMax(viscmax);
    amrex::ParallelDescriptor::ReduceRealMax(elasticmax);
    amrex::ParallelDescriptor::ReduceRealMax(phasefieldmax);

    Set::Scalar adv_dt = cfl_v;
    if (advmax > 0.0) adv_dt = cfl / advmax;
    Set::Scalar visc_dt = viscmax > 0.0 ? 0.5 * cfl / viscmax : cfl_v;
    Set::Scalar elastic_dt = elasticmax > 0.0 ? cfl / elasticmax : cfl_v;
    Set::Scalar phasefield_dt = phasefieldmax > 0.0 ? cfl / phasefieldmax : cfl_v;
    DynamicTimestep_SyncTimeStep(0, std::min({adv_dt, visc_dt, elastic_dt, phasefield_dt}));
    DynamicTimestep_Update();
}


void
LowMach::PrintDiagnostics(Set::Scalar time, int iter)
{
    if (diagnostics_interval <= 0 || iter % diagnostics_interval != 0) return;

    Set::Scalar vmax = 0.0;
    Set::Scalar uxmin = 1.0e300;
    Set::Scalar uxmax = -1.0e300;
    Set::Scalar uymin = 1.0e300;
    Set::Scalar uymax = -1.0e300;
    Set::Scalar divmax = 0.0;
    Set::Scalar div2 = 0.0;
    Set::Scalar ncell = 0.0;
    Set::Scalar divmax_interior = 0.0;
    Set::Scalar div2_interior = 0.0;
    Set::Scalar ncell_interior = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = velocity_mf[lev]->ixType() == amrex::IndexType::TheNodeType() ? mfi.tilebox() : mfi.tilebox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpMin, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                Set::Scalar speed = std::sqrt(u(i,j,k,0) * u(i,j,k,0) + u(i,j,k,1) * u(i,j,k,1));
                return {speed, u(i,j,k,0), u(i,j,k,0), u(i,j,k,1), u(i,j,k,1)};
            });
            ReduceTuple hv = reduce_data.value();
            vmax = std::max(vmax, amrex::get<0>(hv));
            uxmin = std::min(uxmin, amrex::get<1>(hv));
            uxmax = std::max(uxmax, amrex::get<2>(hv));
            uymin = std::min(uymin, amrex::get<3>(hv));
            uymax = std::max(uymax, amrex::get<4>(hv));
        }

        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum,
                              amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar,
                              Set::Scalar, Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                auto sten = Numeric::GetStencil(i, j, k, domain);
                Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
                Set::Scalar div = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) div += grad_u(d,d);
                bool interior = i > domain.smallEnd(0) + 2 && i < domain.bigEnd(0) - 2 &&
                                j > domain.smallEnd(1) + 2 && j < domain.bigEnd(1) - 2;
                Set::Scalar abs_div = Util::Abs(div);
                return {abs_div, div * div, 1.0,
                        interior ? abs_div : 0.0,
                        interior ? div * div : 0.0,
                        interior ? 1.0 : 0.0};
            });
            ReduceTuple hv = reduce_data.value();
            divmax = std::max(divmax, amrex::get<0>(hv));
            div2 += amrex::get<1>(hv);
            ncell += amrex::get<2>(hv);
            divmax_interior = std::max(divmax_interior, amrex::get<3>(hv));
            div2_interior += amrex::get<4>(hv);
            ncell_interior += amrex::get<5>(hv);
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(vmax);
    amrex::ParallelDescriptor::ReduceRealMin(uxmin);
    amrex::ParallelDescriptor::ReduceRealMax(uxmax);
    amrex::ParallelDescriptor::ReduceRealMin(uymin);
    amrex::ParallelDescriptor::ReduceRealMax(uymax);
    amrex::ParallelDescriptor::ReduceRealMax(divmax);
    amrex::ParallelDescriptor::ReduceRealSum(div2);
    amrex::ParallelDescriptor::ReduceRealSum(ncell);
    amrex::ParallelDescriptor::ReduceRealMax(divmax_interior);
    amrex::ParallelDescriptor::ReduceRealSum(div2_interior);
    amrex::ParallelDescriptor::ReduceRealSum(ncell_interior);
    Set::Scalar divrms = ncell > 0.0 ? std::sqrt(div2 / ncell) : 0.0;
    Set::Scalar divrms_interior = ncell_interior > 0.0 ? std::sqrt(div2_interior / ncell_interior) : 0.0;
    if (amrex::ParallelDescriptor::IOProcessor())
        amrex::Print() << "LowMach diagnostics step " << iter
                       << " time " << time
                       << " vmax " << vmax
                       << " uxmin " << uxmin
                       << " uxmax " << uxmax
                       << " uymin " << uymin
                       << " uymax " << uymax
                       << " divmax " << divmax
                       << " divrms " << divrms
                       << " divmax_interior " << divmax_interior
                       << " divrms_interior " << divrms_interior << "\n";
}

void
LowMach::TimeStepComplete(Set::Scalar time, int iter)
{
    ProjectVelocity(time + dt[0], dt[0]);
    PrintDiagnostics(time, iter);
}

void
LowMach::PreparePlotFile(Set::Scalar /*time*/, const amrex::Vector<int>& /*iter*/)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        UpdateThermodynamics(lev, *temperature_mf[lev], *mass_fraction_mf[lev]);
        UpdateSolidStress(lev, *velocity_mf[lev], *temperature_mf[lev], *eta_mf[lev], *xi_mf[lev], diagnostics_extended_fields);
        UpdateDerivedDiagnostics(lev, *velocity_mf[lev], *temperature_mf[lev]);
    }
}

void
LowMach::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = std::sqrt(DX[0] * DX[0] + DX[1] * DX[1]);

    for (amrex::MFIter mfi(*temperature_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tag = tags.array(mfi);
        Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
        const Set::Scalar vcrit = velocity_refinement_criterion;
        const Set::Scalar pcrit = pressure_refinement_criterion;
        const Set::Scalar Tcrit = temperature_refinement_criterion;
        const Set::Scalar etacrit = eta_refinement_criterion;
        amrex::Box domain = geom[lev].Domain();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Vector grad_p = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
            Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX, sten);
            const Set::Scalar eta_val = eta(i,j,k);
            const bool eta_interface = eta_val > etacrit && eta_val < 1.0 - etacrit;
            if (grad_u.norm() * dr > vcrit ||
                grad_p.lpNorm<2>() * dr > pcrit ||
                grad_T.lpNorm<2>() * dr > Tcrit ||
                grad_eta.lpNorm<2>() * dr * 2.0 > etacrit ||
                eta_interface)
                tag(i,j,k) = amrex::TagBox::SET;
        });
    }
}

}
