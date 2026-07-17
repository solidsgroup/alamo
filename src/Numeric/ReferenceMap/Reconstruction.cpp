#include "Numeric/ReferenceMap/Reconstruction.H"

#include <cmath>

#include "AMReX_Gpu.H"
#include "AMReX_MFIter.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Numeric
{
namespace ReferenceMap
{
void
Reconstruction::Parse(Reconstruction& value, IO::ParmParse& pp)
{
    pp.query_default("eta_core", value.eta_core, 0.5);
    pp.query_default("eta_extension", value.eta_extension, 1.0e-3);
    pp.query_default("extrapolation_sweeps", value.extrapolation_sweeps, 4);
    pp.query_default("reconstruction_alpha", value.reconstruction_alpha, 1.0);
    pp.query_default("reconstruction_power", value.reconstruction_power, 1.0);
    pp.query_default("affine_tolerance", value.affine_tolerance, 1.0e-7);
    pp.query_default("smoothing_sweeps", value.smoothing_sweeps, 2);
    pp.query_default("smoothing_alpha", value.smoothing_alpha, 1.0);
    pp.query_default("smoothing_power", value.smoothing_power, 2.0);
    pp.query_default("stress_smoothing_sweeps", value.stress_smoothing_sweeps, 0);
    pp.query_default("stress_smoothing_alpha", value.stress_smoothing_alpha, 0.1);
}

void
Reconstruction::SmoothForStress(const amrex::Geometry& geom,
                                const amrex::MultiFab& eta_mf,
                                amrex::MultiFab& xi_mf) const
{
    const int sweeps = Util::Max(stress_smoothing_sweeps, 0);
    if (sweeps == 0 || stress_smoothing_alpha == 0.0) return;

    const Set::Scalar stress_eta_min = Util::Clamp(eta_extension, 0.0, 1.0);
    const Set::Scalar core_eta = Util::Max(Util::Clamp(eta_core, 0.0, 1.0), 1.0e-12);
    const Set::Scalar power = Util::Max(smoothing_power, 0.0);
    const int xi_ngrow = xi_mf.nGrow();
    const amrex::Box domain = geom.Domain();
    const amrex::Dim3 lo = amrex::lbound(domain);
    const amrex::Dim3 hi = amrex::ubound(domain);

    for (int sweep = 0; sweep < sweeps; ++sweep)
    {
        xi_mf.FillBoundary(geom.periodicity());

        amrex::MultiFab xi_old(
            xi_mf.boxArray(), xi_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        xi_old.FillBoundary(geom.periodicity());

        for (amrex::MFIter mfi(xi_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const amrex::Array4<const Set::Scalar> eta = eta_mf.const_array(mfi);
            const amrex::Array4<const Set::Scalar> xi_in = xi_old.const_array(mfi);
            const amrex::Array4<Set::Scalar> xi = xi_mf.array(mfi);
            const amrex::Geometry geometry = geom;
            const Set::Scalar alpha = stress_smoothing_alpha;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
                if (eta_val <= stress_eta_min) return;

                const Set::Scalar grade = 1.0 - Util::SmootherStep(eta_val / core_eta);
                if (grade == 0.0) return;
                const Set::Scalar local_alpha = alpha * std::pow(grade, power);
                if (local_alpha == 0.0) return;

                const Set::Vector pos =
                    Set::Position(i, j, k, geometry, amrex::IndexType::TheCellType());
                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    const Set::Scalar q = xi_in(i,j,k,n) - pos(n);
                    Set::Scalar lap = 0.0;

                    if (i > lo.x)
                    {
                        const Set::Vector npos = Set::Position(
                            i-1, j, k, geometry, amrex::IndexType::TheCellType());
                        lap += xi_in(i-1,j,k,n) - npos(n) - q;
                    }
                    if (i < hi.x)
                    {
                        const Set::Vector npos = Set::Position(
                            i+1, j, k, geometry, amrex::IndexType::TheCellType());
                        lap += xi_in(i+1,j,k,n) - npos(n) - q;
                    }
#if AMREX_SPACEDIM >= 2
                    if (j > lo.y)
                    {
                        const Set::Vector npos = Set::Position(
                            i, j-1, k, geometry, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j-1,k,n) - npos(n) - q;
                    }
                    if (j < hi.y)
                    {
                        const Set::Vector npos = Set::Position(
                            i, j+1, k, geometry, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j+1,k,n) - npos(n) - q;
                    }
#endif
#if AMREX_SPACEDIM == 3
                    if (k > lo.z)
                    {
                        const Set::Vector npos = Set::Position(
                            i, j, k-1, geometry, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j,k-1,n) - npos(n) - q;
                    }
                    if (k < hi.z)
                    {
                        const Set::Vector npos = Set::Position(
                            i, j, k+1, geometry, amrex::IndexType::TheCellType());
                        lap += xi_in(i,j,k+1,n) - npos(n) - q;
                    }
#endif
                    xi(i,j,k,n) = pos(n) + q + local_alpha * lap;
                }
            });
        }
    }

    xi_mf.FillBoundary(geom.periodicity());
}

void
Reconstruction::operator()(const amrex::Geometry& geom,
                           const amrex::MultiFab& eta_mf,
                           amrex::MultiFab& xi_mf,
                           BC::BC<Set::Scalar>& xi_bc,
                           const Set::Scalar time) const
{
    const Set::Scalar core_eta =
        Util::Max(Util::Clamp(eta_core, 0.0, 1.0), 1.0e-12);
    const Set::Scalar repair_alpha = Util::Clamp(reconstruction_alpha, 0.0, 1.0);
    const Set::Scalar repair_power = Util::Max(reconstruction_power, 0.0);
    const Set::Scalar curvature_tolerance = Util::Max(affine_tolerance, 1.0e-12);
    const Set::Scalar smooth_alpha = Util::Clamp(smoothing_alpha, 0.0, 1.0);
    const Set::Scalar smooth_power = Util::Max(smoothing_power, 0.0);
    const int extrap_sweeps = Util::Max(extrapolation_sweeps, 0);
    const int smooth_sweeps = Util::Max(smoothing_sweeps, 0);
    const int xi_ngrow = xi_mf.nGrow();
    const amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM> DX = geom.CellSizeArray();

    amrex::MultiFab eta_work(
        eta_mf.boxArray(), eta_mf.DistributionMap(), 1, eta_mf.nGrow());
    amrex::MultiFab::Copy(eta_work, eta_mf, 0, 0, 1, eta_mf.nGrow());
    eta_work.FillBoundary(geom.periodicity());

    auto fill_xi_boundary = [&]()
    {
        xi_bc.define(geom);
        xi_bc.FillBoundary(xi_mf, 0, AMREX_SPACEDIM, time, 0);
        xi_mf.FillBoundary(geom.periodicity());
    };

    fill_xi_boundary();

    for (int sweep = 0; sweep < extrap_sweeps; ++sweep)
    {
        amrex::MultiFab xi_old(
            xi_mf.boxArray(), xi_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        xi_old.FillBoundary(geom.periodicity());

        for (amrex::MFIter mfi(xi_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const amrex::Array4<const Set::Scalar> eta = eta_work.const_array(mfi);
            const amrex::Array4<const Set::Scalar> xi_in = xi_old.const_array(mfi);
            const amrex::Array4<Set::Scalar> xi = xi_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
                const Set::Scalar grade = 1.0 - Util::SmootherStep(eta_val / core_eta);
                if (grade == 0.0) return;
                const Set::Scalar repair = repair_alpha * std::pow(grade, repair_power);
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
                            const Set::Scalar eta_near =
                                Util::Clamp(eta(ii,jj,kk), 0.0, 1.0);
                            if (eta_near <= eta_val) continue;

                            const int iii = i + (d == 0 ? 2 * sign : 0);
                            const int jjj = j + (d == 1 ? 2 * sign : 0);
                            const int kkk = k + (d == 2 ? 2 * sign : 0);
                            const Set::Scalar eta_far =
                                Util::Clamp(eta(iii,jjj,kkk), 0.0, 1.0);
                            const Set::Scalar nearest =
                                xi_in(ii,jj,kk,n) - (n == d ? sign * DX[d] : 0.0);
                            Set::Scalar candidate = nearest;
                            if (eta_far > eta_near)
                            {
                                const Set::Scalar inward_slope =
                                    xi_in(ii,jj,kk,n) - xi_in(iii,jjj,kkk,n);
                                const Set::Scalar linear = xi_in(ii,jj,kk,n) + inward_slope;
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
                                    const Set::Scalar affine_weight =
                                        1.0 - Util::SmootherStep(
                                            relative_curvature / curvature_tolerance);
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
                        xi(i,j,k,n) =
                            (1.0 - repair) * xi_in(i,j,k,n) + repair * reconstructed;
                    }
                }
            });
        }
        amrex::Gpu::streamSynchronize();
        fill_xi_boundary();
    }

    for (int sweep = 0; sweep < smooth_sweeps; ++sweep)
    {
        amrex::MultiFab xi_old(
            xi_mf.boxArray(), xi_mf.DistributionMap(), AMREX_SPACEDIM, xi_ngrow);
        amrex::MultiFab::Copy(xi_old, xi_mf, 0, 0, AMREX_SPACEDIM, xi_ngrow);
        xi_old.FillBoundary(geom.periodicity());

        for (amrex::MFIter mfi(xi_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const amrex::Array4<const Set::Scalar> eta = eta_work.const_array(mfi);
            const amrex::Array4<const Set::Scalar> xi_in = xi_old.const_array(mfi);
            const amrex::Array4<Set::Scalar> xi = xi_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
                const Set::Scalar grade = 1.0 - Util::SmootherStep(eta_val / core_eta);
                if (grade == 0.0) return;
                const Set::Scalar relax = smooth_alpha * std::pow(grade, smooth_power);
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
                    xi(i,j,k,n) =
                        (1.0 - relax) * xi_in(i,j,k,n) + relax * sum / count;
                }
            });
        }
        amrex::Gpu::streamSynchronize();
        fill_xi_boundary();
    }

    fill_xi_boundary();
}
}
}
