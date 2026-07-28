#include "Diffusion.H"

#include <limits>

#include <AMReX_FillPatchUtil.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>

#include "BC/Constant.H"

namespace
{
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
TransverseGradient(const amrex::Array4<const Set::Scalar>& state,
                    int i, int j, int k, int n,
                    int face_direction, int derivative_direction,
                    const amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM>& dxinv)
{
    const int ilo = i - (face_direction == 0);
    const int jlo = j - (face_direction == 1);
    const int klo = k - (face_direction == 2);
    const int di = derivative_direction == 0;
    const int dj = derivative_direction == 1;
    const int dk = derivative_direction == 2;
    return 0.25 * dxinv[derivative_direction] *
        (state(i + di, j + dj, k + dk, n) -
         state(i - di, j - dj, k - dk, n) +
         state(ilo + di, jlo + dj, klo + dk, n) -
         state(ilo - di, jlo - dj, klo - dk, n));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
CrossFlux(const amrex::Array4<const Set::Scalar>& state,
          const amrex::GpuArray<amrex::Array4<const Set::Scalar>,
                                AMREX_SPACEDIM>& tensor,
          int i, int j, int k, int n, int face_direction,
          const amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM>& dxinv,
          Set::Scalar scale)
{
    Set::Scalar flux = 0.0;
    for (int derivative_direction = 0;
         derivative_direction < AMREX_SPACEDIM; ++derivative_direction)
        if (derivative_direction != face_direction)
            flux -= scale * tensor[face_direction](
                i,j,k,derivative_direction) *
                TransverseGradient(state, i, j, k, n, face_direction,
                                    derivative_direction, dxinv);
    return flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar
TensorFlux(const amrex::Array4<const Set::Scalar>& state,
           const amrex::GpuArray<amrex::Array4<const Set::Scalar>,
                                 AMREX_SPACEDIM>& tensor,
           int i, int j, int k, int n, int face_direction,
           const amrex::GpuArray<Set::Scalar, AMREX_SPACEDIM>& dxinv,
           Set::Scalar scale)
{
    const int ilo = i - (face_direction == 0);
    const int jlo = j - (face_direction == 1);
    const int klo = k - (face_direction == 2);
    return -scale * tensor[face_direction](
               i,j,k,face_direction) * dxinv[face_direction] *
               (state(i,j,k,n) - state(ilo,jlo,klo,n)) +
           CrossFlux(state, tensor, i, j, k, n, face_direction,
                     dxinv, scale);
}

class MLTensorDiffusion : public amrex::MLABecLaplacian
{
public:
    MLTensorDiffusion(
        const amrex::Vector<amrex::Geometry>& geometry,
        const amrex::Vector<amrex::BoxArray>& grids,
        const amrex::Vector<amrex::DistributionMapping>& distribution_mapping,
        const amrex::LPInfo& info, int ncomp)
        : amrex::MLABecLaplacian(
            geometry, grids, distribution_mapping, info, {}, ncomp)
    {
        tensor_coefficients.resize(m_num_amr_levels);
        for (int alev = 0; alev < m_num_amr_levels; ++alev)
        {
            tensor_coefficients[alev].resize(m_num_mg_levels[alev]);
            for (int mglev = 0; mglev < m_num_mg_levels[alev]; ++mglev)
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                {
                    amrex::BoxArray faces = m_grids[alev][mglev];
                    faces.surroundingNodes(d);
                    tensor_coefficients[alev][mglev][d].define(
                        faces, m_dmap[alev][mglev], AMREX_SPACEDIM, 0);
                }
        }
    }

    bool isCrossStencil() const override { return false; }
    bool isTensorOp() const override { return true; }

    void SetTensorCoefficients(
        int alev,
        const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& coefficients)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            amrex::MultiFab::Copy(
                tensor_coefficients[alev][0][d], *coefficients[d],
                0, 0, AMREX_SPACEDIM, 0);
    }

    void prepareForSolve() override
    {
        amrex::MLABecLaplacian::prepareForSolve();
        for (int alev = m_num_amr_levels - 1; alev > 0; --alev)
        {
            AverageDown(alev);
            amrex::average_down_faces(
                amrex::GetArrOfConstPtrs(tensor_coefficients[alev].back()),
                amrex::GetArrOfPtrs(tensor_coefficients[alev - 1].front()),
                amrex::IntVect(mg_coarsen_ratio), m_geom[alev - 1][0]);
        }
        AverageDown(0);
    }

    void Fapply(int alev, int mglev, amrex::MultiFab& out,
                const amrex::MultiFab& in) const override
    {
        amrex::MLABecLaplacian::Fapply(alev, mglev, out, in);
        AddCrossOperator(alev, mglev, out, in, 1.0);
    }

    void Fsmooth(int alev, int mglev, amrex::MultiFab& state,
                 const amrex::MultiFab& rhs, int redblack) const override
    {
        amrex::ignore_unused(redblack);
        amrex::MultiFab applied(
            state.boxArray(), state.DistributionMap(), state.nComp(), 0);
        Fapply(alev, mglev, applied, state);

        const amrex::MultiFab& mass = m_a_coeffs[alev][mglev];
        const auto dxinv = m_geom[alev][mglev].InvCellSizeArray();
        const Set::Scalar alpha = m_a_scalar;
        const Set::Scalar beta = m_b_scalar;
        constexpr Set::Scalar omega = 0.6;
        for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto x = state.array(mfi);
            const auto Ax = applied.const_array(mfi);
            const auto f = rhs.const_array(mfi);
            const auto a = mass.const_array(mfi);
            const amrex::GpuArray<amrex::Array4<const Set::Scalar>,
                                  AMREX_SPACEDIM> diagonal{{AMREX_D_DECL(
                m_b_coeffs[alev][mglev][0].const_array(mfi),
                m_b_coeffs[alev][mglev][1].const_array(mfi),
                m_b_coeffs[alev][mglev][2].const_array(mfi))}};
            const int ncomp = getNComp();
            amrex::ParallelFor(
                bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                {
                    Set::Scalar diag = alpha * a(i,j,k);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        const int di = d == 0;
                        const int dj = d == 1;
                        const int dk = d == 2;
                        diag += beta * dxinv[d] * dxinv[d] *
                            (diagonal[d](i,j,k,n) +
                             diagonal[d](i + di, j + dj, k + dk, n));
                    }
                    x(i,j,k,n) += omega *
                        (f(i,j,k,n) - Ax(i,j,k,n)) / diag;
                });
        }
    }

    void FFlux(
        int alev, const amrex::MFIter& mfi,
        const amrex::Array<amrex::FArrayBox*, AMREX_SPACEDIM>& flux,
        const amrex::FArrayBox& state, Location location,
        int face_only = 0) const override
    {
        amrex::ignore_unused(location);
        const int mglev = 0;
        const amrex::Box& bx = mfi.tilebox();
        const auto x = state.const_array();
        const auto dxinv = m_geom[alev][mglev].InvCellSizeArray();
        const Set::Scalar beta = m_b_scalar;
        const int ncomp = getNComp();
        const amrex::GpuArray<amrex::Array4<const Set::Scalar>,
                              AMREX_SPACEDIM> tensor{{AMREX_D_DECL(
            tensor_coefficients[alev][mglev][0].const_array(mfi),
            tensor_coefficients[alev][mglev][1].const_array(mfi),
            tensor_coefficients[alev][mglev][2].const_array(mfi))}};

        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            const auto q = flux[d]->array();
            if (face_only)
            {
                const amrex::Box lo = amrex::bdryLo(bx, d);
                const int length = bx.length(d);
                amrex::ParallelFor(
                    lo, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                    {
                        q(i,j,k,n) = TensorFlux(
                            x, tensor, i, j, k, n, d, dxinv, beta);
                        q(i + (d == 0) * length,
                          j + (d == 1) * length,
                          k + (d == 2) * length,n) = TensorFlux(
                            x, tensor,
                            i + (d == 0) * length,
                            j + (d == 1) * length,
                            k + (d == 2) * length,
                            n, d, dxinv, beta);
                    });
            }
            else
            {
                const amrex::Box faces = amrex::surroundingNodes(bx, d);
                amrex::ParallelFor(
                    faces, ncomp,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                    {
                        q(i,j,k,n) = TensorFlux(
                            x, tensor, i, j, k, n, d, dxinv, beta);
                    });
            }
        }
    }

private:
    void AverageDown(int alev)
    {
        for (int mglev = 1; mglev < m_num_mg_levels[alev]; ++mglev)
        {
            const amrex::IntVect ratio = alev > 0 ?
                amrex::IntVect(mg_coarsen_ratio) :
                mg_coarsen_ratio_vec[mglev - 1];
            amrex::average_down_faces(
                amrex::GetArrOfConstPtrs(
                    tensor_coefficients[alev][mglev - 1]),
                amrex::GetArrOfPtrs(tensor_coefficients[alev][mglev]),
                ratio, 0);
        }
    }

    void AddCrossOperator(int alev, int mglev, amrex::MultiFab& out,
                          const amrex::MultiFab& state,
                          Set::Scalar scale) const
    {
        const auto dxinv = m_geom[alev][mglev].InvCellSizeArray();
        const Set::Scalar beta = m_b_scalar;
        const int ncomp = getNComp();
        for (amrex::MFIter mfi(out, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto y = out.array(mfi);
            const auto x = state.const_array(mfi);
            const amrex::GpuArray<amrex::Array4<const Set::Scalar>,
                                  AMREX_SPACEDIM> tensor{{AMREX_D_DECL(
                tensor_coefficients[alev][mglev][0].const_array(mfi),
                tensor_coefficients[alev][mglev][1].const_array(mfi),
                tensor_coefficients[alev][mglev][2].const_array(mfi))}};
            amrex::ParallelFor(
                bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                {
                    Set::Scalar divergence = 0.0;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        const int di = d == 0;
                        const int dj = d == 1;
                        const int dk = d == 2;
                        divergence += dxinv[d] *
                            (CrossFlux(x, tensor, i + di, j + dj, k + dk,
                                       n, d, dxinv, beta) -
                             CrossFlux(x, tensor, i, j, k, n, d,
                                       dxinv, beta));
                    }
                    y(i,j,k,n) += scale * divergence;
                });
        }
    }

    amrex::Vector<amrex::Vector<
        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>>>
        tensor_coefficients;
};
}

namespace Operator
{
void
Diffusion::Parse(Diffusion& value, IO::ParmParse& pp)
{
    pp.query_default("tol_rel", value.tolerance_relative, 1.0e-11);
    pp.query_default("tol_abs", value.tolerance_absolute, 1.0e-12);
    pp.query_default("verbose", value.verbose, 0);
    pp.query_default("max_order", value.max_order, 2);
}

Diffusion::System&
Diffusion::GetSystem(int ncomp)
{
    auto system = systems.find(ncomp);
    Util::Assert(INFO, TEST(system != systems.end()));
    return *system->second;
}

void
Diffusion::SetLayout(
    const amrex::Vector<amrex::Geometry>& a_geometry,
    const amrex::Vector<amrex::IntVect>& a_refinement_ratio,
    const Set::Field<Set::Scalar>& layout,
    int number_of_levels, int number_of_components)
{
    Util::Assert(INFO, TEST(number_of_levels > 0));
    Util::Assert(INFO, TEST(number_of_components > 0));
    Util::Assert(INFO, TEST(number_of_levels <= static_cast<int>(a_geometry.size())));
    Util::Assert(INFO, TEST(number_of_levels <= static_cast<int>(layout.size())));
    Util::Assert(INFO, TEST(number_of_levels == 1 ||
                            number_of_levels - 1 <=
                                static_cast<int>(a_refinement_ratio.size())));

    nlevels = number_of_levels;
    geometry.assign(a_geometry.begin(), a_geometry.begin() + nlevels);
    refinement_ratio.assign(
        a_refinement_ratio.begin(), a_refinement_ratio.begin() + nlevels - 1);

    auto& system_pointer = systems[number_of_components];
    bool layout_changed = !system_pointer ||
        static_cast<int>(system_pointer->state.size()) != nlevels;
    if (!layout_changed)
    {
        for (int lev = 0; lev < nlevels; ++lev)
        {
            if (system_pointer->state[lev]->boxArray() != layout[lev]->boxArray() ||
                system_pointer->state[lev]->DistributionMap() !=
                    layout[lev]->DistributionMap())
            {
                layout_changed = true;
                break;
            }
        }
    }
    if (!layout_changed) return;

    grids.resize(nlevels);
    distribution_mapping.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        grids[lev] = layout[lev]->boxArray();
        distribution_mapping[lev] = layout[lev]->DistributionMap();
    }

    system_pointer = std::make_unique<System>();
    System& system = *system_pointer;
    system.state.Define(nlevels, grids, distribution_mapping,
                        number_of_components, 1);
    system.source.Define(nlevels, grids, distribution_mapping,
                        number_of_components, 0);
    system.rhs.Define(nlevels, grids, distribution_mapping,
                        number_of_components, 0);
    system.mass.Define(nlevels, grids, distribution_mapping, 1, 0);
    system.mobility.Define(nlevels, grids, distribution_mapping, 1, 1);
    system.face_mobility.resize(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            amrex::BoxArray face_grids = grids[lev];
            face_grids.surroundingNodes(d);
            system.face_mobility[lev][d].define(
                face_grids, distribution_mapping[lev], 1, 0);
        }
    }
}

amrex::MultiFab&
Diffusion::State(int lev, int ncomp)
{
    return *GetSystem(ncomp).state[lev];
}

amrex::MultiFab&
Diffusion::Source(int lev, int ncomp)
{
    return *GetSystem(ncomp).source[lev];
}

amrex::MultiFab&
Diffusion::Mass(int lev, int ncomp)
{
    return *GetSystem(ncomp).mass[lev];
}

amrex::MultiFab&
Diffusion::Mobility(int lev, int ncomp)
{
    return *GetSystem(ncomp).mobility[lev];
}

amrex::MultiFab&
Diffusion::TensorMobility(int lev, int ncomp)
{
    System& system = GetSystem(ncomp);
    if (system.tensor_mobility.size() == 0)
    {
        system.tensor_mobility.Define(
            nlevels, grids, distribution_mapping,
            AMREX_SPACEDIM * AMREX_SPACEDIM, 1);
        system.face_tensor_mobility.resize(nlevels);
        for (int level = 0; level < nlevels; ++level)
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                amrex::BoxArray face_grids = grids[level];
                face_grids.surroundingNodes(d);
                system.face_tensor_mobility[level][d].define(
                    face_grids, distribution_mapping[level],
                    AMREX_SPACEDIM, 0);
            }
    }
    return *system.tensor_mobility[lev];
}

void
Diffusion::Solve(Set::Scalar time, Set::Scalar dt,
                const amrex::BCRec& boundary_condition, int ncomp,
                bool use_tensor_mobility, bool include_source)
{
    Solve(time, dt, amrex::Vector<amrex::BCRec>(ncomp, boundary_condition),
        ncomp, use_tensor_mobility, include_source);
}

void
Diffusion::Solve(Set::Scalar time, Set::Scalar dt,
                const amrex::Vector<amrex::BCRec>& boundary_conditions,
                int ncomp, bool use_tensor_mobility, bool include_source)
{
    BL_PROFILE("Operator::Diffusion::Solve");
    Util::Assert(INFO,
        TEST(static_cast<int>(boundary_conditions.size()) == ncomp));
    System& system = GetSystem(ncomp);
    if (use_tensor_mobility)
        Util::Assert(INFO, TEST(system.tensor_mobility.size() > 0));
    BC::Constant::ZeroNeumann coefficient_bc(1);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        coefficient_bc.define(geometry[lev]);
        if (lev == 0)
        {
            coefficient_bc.FillBoundary(
                *system.mobility[lev], 0, 1, time, 0);
            system.mobility[lev]->FillBoundary(geometry[lev].periodicity());
        }
        else
        {
            amrex::Vector<amrex::MultiFab*> coarse{
                system.mobility[lev - 1].get()};
            amrex::Vector<amrex::MultiFab*> fine{
                system.mobility[lev].get()};
            amrex::Vector<amrex::Real> coarse_time{time};
            amrex::Vector<amrex::Real> fine_time{time};
            amrex::Vector<amrex::BCRec> bcs(1, coefficient_bc.GetBCRec());
            amrex::FillPatchTwoLevels(
                *system.mobility[lev], time,
                coarse, coarse_time, fine, fine_time,
                0, 0, 1,
                geometry[lev - 1], geometry[lev],
                coefficient_bc, 0, coefficient_bc, 0,
                refinement_ratio[lev - 1], &amrex::cell_cons_interp, bcs, 0);
        }
    }

    if (use_tensor_mobility)
    {
        constexpr int tensor_components =
            AMREX_SPACEDIM * AMREX_SPACEDIM;
        BC::Constant::ZeroNeumann tensor_bc(tensor_components);
        for (int lev = 0; lev < nlevels; ++lev)
        {
            tensor_bc.define(geometry[lev]);
            if (lev == 0)
            {
                tensor_bc.FillBoundary(
                    *system.tensor_mobility[lev], 0,
                    tensor_components, time, 0);
                system.tensor_mobility[lev]->FillBoundary(
                    geometry[lev].periodicity());
            }
            else
            {
                amrex::Vector<amrex::MultiFab*> coarse{
                    system.tensor_mobility[lev - 1].get()};
                amrex::Vector<amrex::MultiFab*> fine{
                    system.tensor_mobility[lev].get()};
                amrex::Vector<amrex::Real> coarse_time{time};
                amrex::Vector<amrex::Real> fine_time{time};
                amrex::Vector<amrex::BCRec> bcs(
                    tensor_components, tensor_bc.GetBCRec());
                amrex::FillPatchTwoLevels(
                    *system.tensor_mobility[lev], time,
                    coarse, coarse_time, fine, fine_time,
                    0, 0, tensor_components,
                    geometry[lev - 1], geometry[lev],
                    tensor_bc, 0, tensor_bc, 0,
                    refinement_ratio[lev - 1],
                    &amrex::cell_cons_interp, bcs, 0);
            }
        }
    }

    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>>
        face_mobility_pointer(nlevels);
    amrex::Vector<amrex::Array<amrex::MultiFab const*, AMREX_SPACEDIM>>
        face_tensor_mobility_pointer(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> face_pointer;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            face_pointer[d] = &system.face_mobility[lev][d];
            face_mobility_pointer[lev][d] = &system.face_mobility[lev][d];
            if (use_tensor_mobility)
                face_tensor_mobility_pointer[lev][d] =
                    &system.face_tensor_mobility[lev][d];
        }
        amrex::average_cellcenter_to_face(
            face_pointer, *system.mobility[lev], geometry[lev], 1, true, 0);
        if (use_tensor_mobility)
        {
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                for (amrex::MFIter mfi(system.face_tensor_mobility[lev][d],
                                        amrex::TilingIfNotGPU());
                    mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.tilebox();
                    const auto scalar_cell =
                        system.mobility[lev]->const_array(mfi);
                    const auto tensor_cell =
                        system.tensor_mobility[lev]->const_array(mfi);
                    const auto tensor_face =
                        system.face_tensor_mobility[lev][d].array(mfi);
                    const auto diagonal_face =
                        system.face_mobility[lev][d].array(mfi);
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            const int ilo = i - (d == 0);
                            const int jlo = j - (d == 1);
                            const int klo = k - (d == 2);
                            Set::Matrix lo = Set::Matrix::Zero();
                            Set::Matrix hi = Set::Matrix::Zero();
                            for (int e = 0; e < AMREX_SPACEDIM; ++e)
                                for (int f = 0; f < AMREX_SPACEDIM; ++f)
                                {
                                    lo(e,f) = tensor_cell(
                                        ilo,jlo,klo,
                                        e * AMREX_SPACEDIM + f);
                                    hi(e,f) = tensor_cell(
                                        i,j,k,e * AMREX_SPACEDIM + f);
                                }
                            lo.diagonal().array() +=
                                scalar_cell(ilo,jlo,klo);
                            hi.diagonal().array() += scalar_cell(i,j,k);
                            const Set::Matrix face =
                                2.0 * (lo.inverse() + hi.inverse()).inverse();
                            for (int e = 0; e < AMREX_SPACEDIM; ++e)
                                tensor_face(i,j,k,e) = face(d,e);
                            diagonal_face(i,j,k) = face(d,d);
                        });
                }
            }
        }

        for (amrex::MFIter mfi(*system.rhs[lev], amrex::TilingIfNotGPU());
            mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const auto state = system.state[lev]->const_array(mfi);
            const auto source = system.source[lev]->const_array(mfi);
            const auto mass = system.mass[lev]->const_array(mfi);
            const auto rhs = system.rhs[lev]->array(mfi);
            amrex::ParallelFor(
                bx, ncomp,
                [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                {
                    rhs(i,j,k,n) = mass(i,j,k) * state(i,j,k,n) +
                        (include_source ? dt * source(i,j,k,n) : 0.0);
                });
        }
    }

    amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
        boundary_lo(ncomp);
    amrex::Vector<amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>>
        boundary_hi(ncomp);
    for (int n = 0; n < ncomp; ++n)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            if (geometry[0].isPeriodic(d))
            {
                boundary_lo[n][d] = amrex::LinOpBCType::Periodic;
                boundary_hi[n][d] = amrex::LinOpBCType::Periodic;
            }
            else
            {
                boundary_lo[n][d] =
                    BC::BCUtil::IsDirichlet(boundary_conditions[n].lo(d)) ?
                    amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
                boundary_hi[n][d] =
                    BC::BCUtil::IsDirichlet(boundary_conditions[n].hi(d)) ?
                    amrex::LinOpBCType::Dirichlet : amrex::LinOpBCType::Neumann;
            }
        }
    }

    int semicoarsening_direction = 0;
    Set::Scalar mobility_max = 0.0;
    Set::Scalar mobility_min = std::numeric_limits<Set::Scalar>::max();
    if (use_tensor_mobility)
    {
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            Set::Scalar directional_max = 0.0;
            for (int lev = 0; lev < nlevels; ++lev)
                directional_max = Util::Max(
                    directional_max,
                    system.face_mobility[lev][d].max(0));
            if (directional_max > mobility_max)
            {
                mobility_max = directional_max;
                semicoarsening_direction = d;
            }
            mobility_min = Util::Min(mobility_min, directional_max);
        }
    }
    const bool anisotropic = use_tensor_mobility &&
        mobility_max > 4.0 * mobility_min;
    // AMReX's 2D ABec line smoother currently supports semicoarsening only
    // with the y direction left uncoarsened.
    const bool semicoarsening = anisotropic &&
        semicoarsening_direction == 1;

    amrex::LPInfo info;
    if (anisotropic)
        info.setSemicoarsening(true)
            .setMaxSemicoarseningLevel(100)
            .setSemicoarseningDirection(semicoarsening_direction);
    std::unique_ptr<amrex::MLABecLaplacian> diffusion;
    MLTensorDiffusion* tensor_diffusion = nullptr;
    if (use_tensor_mobility)
    {
        auto tensor = std::make_unique<MLTensorDiffusion>(
            geometry, grids, distribution_mapping, info, ncomp);
        tensor_diffusion = tensor.get();
        diffusion = std::move(tensor);
    }
    else
        diffusion = std::make_unique<amrex::MLABecLaplacian>(
            geometry, grids, distribution_mapping, info,
            amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>{},
            ncomp);

    diffusion->setMaxOrder(max_order);
    diffusion->setDomainBC(boundary_lo, boundary_hi);
    diffusion->setScalars(1.0, dt);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        diffusion->setLevelBC(lev, system.state[lev].get());
        diffusion->setACoeffs(lev, *system.mass[lev]);
        diffusion->setBCoeffs(lev, face_mobility_pointer[lev]);
        if (tensor_diffusion)
            tensor_diffusion->SetTensorCoefficients(
                lev, face_tensor_mobility_pointer[lev]);
    }

    amrex::Vector<amrex::MultiFab*> state_pointer(nlevels);
    amrex::Vector<amrex::MultiFab const*> rhs_pointer(nlevels);
    for (int lev = 0; lev < nlevels; ++lev)
    {
        state_pointer[lev] = system.state[lev].get();
        rhs_pointer[lev] = system.rhs[lev].get();
    }

    amrex::MLMG solver(*diffusion);
    solver.setVerbose(verbose);
    solver.setFinalFillBC(true);
    if (use_tensor_mobility || semicoarsening)
    {
        solver.setPreSmooth(8);
        solver.setPostSmooth(8);
    }
    solver.solve(state_pointer, rhs_pointer,
                tolerance_relative, tolerance_absolute);
}

void
Diffusion::FillBoundary(Set::Field<Set::Scalar>& field,
                        BC::BC<Set::Scalar>& boundary_condition,
                        Set::Scalar time, int ncomp) const
{
    for (int lev = 0; lev < nlevels; ++lev)
    {
        boundary_condition.define(geometry[lev]);
        if (lev == 0)
        {
            amrex::Vector<amrex::MultiFab*> source{field[lev].get()};
            amrex::Vector<amrex::Real> source_time{time};
            amrex::FillPatchSingleLevel(
                *field[lev], time, source, source_time,
                0, 0, ncomp, geometry[lev], boundary_condition, 0);
        }
        else
        {
            amrex::Vector<amrex::MultiFab*> coarse{field[lev - 1].get()};
            amrex::Vector<amrex::MultiFab*> fine{field[lev].get()};
            amrex::Vector<amrex::Real> coarse_time{time};
            amrex::Vector<amrex::Real> fine_time{time};
            amrex::Vector<amrex::BCRec> bcs(ncomp);
            for (int n = 0; n < ncomp; ++n)
                bcs[n] = boundary_condition.GetBCRec(n);
            amrex::FillPatchTwoLevels(
                *field[lev], time,
                coarse, coarse_time, fine, fine_time,
                0, 0, ncomp, geometry[lev - 1], geometry[lev],
                boundary_condition, 0, boundary_condition, 0,
                refinement_ratio[lev - 1], &amrex::cell_cons_interp, bcs, 0);
        }
    }
}

void
Diffusion::Synchronize(Set::Field<Set::Scalar>& field, int ncomp) const
{
    for (int lev = nlevels - 2; lev >= 0; --lev)
        amrex::average_down(*field[lev + 1], *field[lev],
                            geometry[lev + 1], geometry[lev],
                            0, ncomp, refinement_ratio[lev]);
}
}
