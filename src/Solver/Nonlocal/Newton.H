#ifndef SOLVER_NONLOCAL_NEWTON
#define SOLVER_NONLOCAL_NEWTON

#include <cmath>
#include "Set/Set.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"

namespace Solver
{
namespace Nonlocal
{
namespace Detail
{
template <typename T>
inline T
FaceModel(const amrex::Array4<const T>& model,
    const int i, const int j, const int k, const int face)
{
    const int ip = i + (face == 0);
    const int jp = j + (face == 1);
    const int kp = k + (face == 2);
    // Use the same arithmetic face model for the residual and tangent. Tensor
    // models do not have a component-independent harmonic average.
    return 0.5 * (model(i, j, k) + model(ip, jp, kp));
}

template <typename T, typename U>
inline Set::Matrix
FaceKinematicVariable(const T& model,
    const U& u,
    const int i, const int j, const int k, const int face,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix gradu = Numeric::FaceGradient(u, i, j, k, face, dx);
    if (model.kinvar == Model::Solid::KinematicVariable::gradu)
        return gradu;
    if (model.kinvar == Model::Solid::KinematicVariable::epsilon)
        return 0.5 * (gradu + gradu.transpose());
    if (model.kinvar == Model::Solid::KinematicVariable::F)
        return gradu + Set::Matrix::Identity();
    return Set::Matrix::Zero();
}

template <typename T, typename U>
inline Set::Matrix
FaceStress(const amrex::Array4<const T>& model,
    const U& u,
    const int i, const int j, const int k, const int face,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    const T face_model = FaceModel(model, i, j, k, face);
    return face_model.DW(
        FaceKinematicVariable(face_model, u, i, j, k, face, dx));
}

template <typename T, typename U>
inline Set::Matrix4<AMREX_SPACEDIM, T::sym>
FaceTangent(const amrex::Array4<const T>& model,
    const U& u,
    const int i, const int j, const int k, const int face,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    const T face_model = FaceModel(model, i, j, k, face);
    return face_model.DDW(
        FaceKinematicVariable(face_model, u, i, j, k, face, dx));
}

template <typename T, typename U>
inline Set::Vector
FaceStressDivergence(const amrex::Array4<const T>& model,
    const U& u,
    const int i, const int j, const int k,
    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Vector ret = Set::Vector::Zero();
    for (int face = 0; face < AMREX_SPACEDIM; ++face)
    {
        const int im = i - (face == 0);
        const int jm = j - (face == 1);
        const int km = k - (face == 2);
        const Set::Matrix hi = FaceStress(model, u, i, j, k, face, dx);
        const Set::Matrix lo = FaceStress(model, u, im, jm, km, face, dx);
        ret += (hi.col(face) - lo.col(face)) / dx[face];
    }
    return ret;
}
} // namespace Detail

template <typename T>
class Newton: public Linear
{
public:
    Newton() {};

    Newton(Operator::Elastic<T::sym>& a_op)
    {
        this->Define(a_op);
    };

    ~Newton()
    {
        if (m_defined) Clear();
    }

    void Define(Operator::Elastic<T::sym>& a_op)
    {
        Linear::Define(a_op);
        m_elastic = dynamic_cast<Operator::Elastic<T::sym> *>(linop);
        //m_bc = &m_elastic->GetBC();
    }
    void Clear()
    {
        Linear::Clear();
        m_elastic = nullptr;
        //m_bc = nullptr;
    }


    void setNRIters(int a_nriters) { m_nriters = a_nriters; }

    void setPsi(Set::Field<Set::Scalar>& a_psi)
    {
        m_psi = &a_psi;
        if (m_conservative_face_flux && !m_reported_psi_face_flux_conflict)
        {
            Util::Warning(INFO, "solver.conservative_face_flux is disabled because a psi field is active");
            m_reported_psi_face_flux_conflict = true;
        }
    }

    [[nodiscard]] bool usesConservativeFaceFlux() const
    {
        return m_conservative_face_flux && m_psi == nullptr;
    }

private:
    void prepareForSolve(const Set::Field<Set::Scalar>& a_u_mf,
        const Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<Set::Scalar>& a_rhs_mf,
        Set::Field<Set::Matrix>& a_dw_mf,
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>>& a_ddw_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(1));
            domain.convert(amrex::IntVect::TheNodeVector());
            const Set::Scalar* dx = linop->Geom(lev).CellSize();
            Set::Vector DX(linop->Geom(lev).CellSize());
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                bx = bx & domain;

                amrex::Array4<const T>           const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, T::sym>>  const& ddw = a_ddw_mf[lev]->array(mfi);

                // Set model internal dw and ddw.
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);

                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                    {
                        dw(i, j, k) = model(i, j, k).DW(gradu);
                        ddw(i, j, k) = model(i, j, k).DDW(gradu);
                    }
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                    {
                        Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                        dw(i, j, k) = model(i, j, k).DW(eps);
                        ddw(i, j, k) = model(i, j, k).DDW(eps);
                    }
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Matrix F = gradu + Set::Matrix::Identity();
                        dw(i, j, k) = model(i, j, k).DW(F);
                        ddw(i, j, k) = model(i, j, k).DDW(F);
                    }
                });
            }

            Util::RealFillBoundary(*a_dw_mf[lev], m_elastic->Geom(lev));
            Util::RealFillBoundary(*a_ddw_mf[lev], m_elastic->Geom(lev));
        }

        m_elastic->SetModel(a_ddw_mf);

        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(1));
            domain.convert(amrex::IntVect::TheNodeVector());
            const Set::Scalar* dx = linop->Geom(lev).CellSize();
            const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                bx = bx & domain;
                amrex::Array4<const Set::Scalar>  const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar>  const& b = a_b_mf[lev]->array(mfi);
                amrex::Array4<const Set::Matrix>  const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>        const& rhs = a_rhs_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    // Do this if on the domain boundary
                    if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                    {
                        Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
                        Set::Vector U(AMREX_D_DECL(u(i, j, k, 0), u(i, j, k, 1), u(i, j, k, 2)));
                        Set::Vector ret = m_elastic->GetBC()(U, gradu, dw(i, j, k), i, j, k, bx);
                        for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k, d) - ret(d);
                    }
                    else
                    {
                        Set::Vector divdw = Numeric::Divergence(dw, i, j, k, dx);
                        for (int p = 0; p < AMREX_SPACEDIM; p++) rhs(i, j, k, p) = b(i, j, k, p) - divdw(p);
                    }
                });
            }
            Util::RealFillBoundary(*a_ddw_mf[lev], m_elastic->Geom(lev));
            Util::RealFillBoundary(*a_rhs_mf[lev], m_elastic->Geom(lev));
        }
    }

    void prepareForSolve(const Set::Field<Set::Vector>& a_u_mf,
        const Set::Field<Set::Vector>& a_b_mf,
        Set::Field<Set::Scalar>& a_rhs_mf,
        Set::Field<Set::Matrix>& a_dw_mf,
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>>& a_ddw_mf,
        Set::Field<T>& a_model_mf)
    {
        const bool conservative_solve = usesConservativeFaceFlux();
        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            amrex::Box row_domain(linop->Geom(lev).growPeriodicDomain(1));
            row_domain.convert(amrex::IntVect::TheNodeVector());
            const Set::Scalar* dx = linop->Geom(lev).CellSize();
            Set::Vector DX(linop->Geom(lev).CellSize());
            const amrex::Dim3 hi = amrex::ubound(domain);

            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                bx = bx & domain;

                amrex::Array4<const T>           const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector> const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, T::sym>>  const& ddw = a_ddw_mf[lev]->array(mfi);

                if (conservative_solve)
                {
                    amrex::Box cbx = mfi.grownnodaltilebox(-1, 1) & row_domain;
                    amrex::LoopConcurrentOnCpu(cbx, [=] (int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, cbx);
                        Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
                        Set::Matrix kinvar;
                        if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                            kinvar = gradu;
                        else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                            kinvar = 0.5 * (gradu + gradu.transpose());
                        else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                            kinvar = gradu + Set::Matrix::Identity();
                        else
                            Util::Abort(INFO, "Invalid model kinematic variable");

                        dw(i, j, k) = model(i, j, k).DW(kinvar);
                        ddw(i, j, k, 0) = model(i, j, k).DDW(kinvar);
                        for (int face = 0; face < AMREX_SPACEDIM; ++face)
                        {
                            const int index[3] = {i, j, k};
                            const int upper[3] = {hi.x, hi.y, hi.z};
                            if (index[face] < upper[face])
                                ddw(i, j, k, face + 1) = Detail::FaceTangent(
                                    model, u, i, j, k, face, dx);
                            else
                                ddw(i, j, k, face + 1) = ddw(i, j, k, 0);
                        }
                    });
                    continue;
                }

                // Set model internal dw and ddw.
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
                    Set::Matrix kinvar;
                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                        kinvar = gradu; // gradu
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        kinvar = 0.5 * (gradu + gradu.transpose()); // epsilon
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                        kinvar = gradu + Set::Matrix::Identity(); // F

                    dw(i, j, k) = model(i, j, k).DW(kinvar);
                    ddw(i, j, k) = model(i, j, k).DDW(kinvar);

                });
            }

            a_dw_mf[lev]->setMultiGhost(true);
            a_dw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
            a_ddw_mf[lev]->setMultiGhost(true);
            a_ddw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
        }

        m_elastic->SetModel(a_ddw_mf);
        if (m_psi)
            for (int i = 0; i <= a_b_mf.finest_level; i++)
                m_elastic->SetPsi(i, *(*m_psi)[i]);


        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            amrex::Box row_domain(linop->Geom(lev).growPeriodicDomain(1));
            row_domain.convert(amrex::IntVect::TheNodeVector());
            const Set::Scalar* dx = linop->Geom(lev).CellSize();
            const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::Array4<const Set::Vector>  const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector>  const& b = a_b_mf[lev]->array(mfi);
                amrex::Array4<const Set::Matrix>  const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>        const& rhs = a_rhs_mf[lev]->array(mfi);
                if (conservative_solve)
                {
                    amrex::Box cbx = mfi.grownnodaltilebox(-1, 1) & row_domain;
                    amrex::Array4<const T> const& model = a_model_mf[lev]->array(mfi);
                    amrex::LoopConcurrentOnCpu(cbx, [=] (int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, cbx);
                        const int index[3] = {i, j, k};
                        const int lower[3] = {lo.x, lo.y, lo.z};
                        const int upper[3] = {hi.x, hi.y, hi.z};
                        bool on_boundary = false;
                        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
                            on_boundary = on_boundary ||
                                index[dir] == lower[dir] ||
                                index[dir] == upper[dir];
                        if (on_boundary)
                        {
                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
                            Set::Vector ret = m_elastic->GetBC()(
                                u(i, j, k), gradu, dw(i, j, k), i, j, k, cbx);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Vector divdw = Detail::FaceStressDivergence(
                                model, u, i, j, k, dx);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - divdw(d);
                        }
                    });
                    continue;
                }

                // This is for if psi is not being used and has not been set
                if (m_psi == nullptr)
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, bx);
                        // Do this if on the domain boundary
                        if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                        {
                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
                            Set::Vector ret = m_elastic->GetBC()(u(i, j, k), gradu, dw(i, j, k), i, j, k, bx);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Vector divdw = Numeric::Divergence(dw, i, j, k, dx, sten);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - divdw(d);
                        }
                    });
                }
                // This is EXACTLY THE SAME AS THE ABOVE, excpet that we are now accounting
                // for psi terms. The reason we do this is because there's no good way (yet) to 
                // initialize patches if psi = nullptr
                else
                {
                    const amrex::Dim3 boxlo = amrex::lbound(bx), boxhi = amrex::ubound(bx);
                    amrex::Array4<const Set::Scalar>  const& psi = (*m_psi)[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, bx);
                        // Do this if on the domain boundary
                        if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                        {
                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, dx, sten);
                            Set::Scalar psiavg = 1.0;
                            if (AMREX_D_TERM(i > boxlo.x && i<boxhi.x, && j>boxlo.y && j<boxhi.y, && k>boxlo.z && k < boxhi.z))
                            {
                                psiavg = Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0);
                            }
                            Set::Vector ret = m_elastic->GetBC()(u(i, j, k), gradu, dw(i, j, k) * psiavg, i, j, k, bx);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Vector divdw = Numeric::Divergence(dw, i, j, k, dx, sten);
                            if (AMREX_D_TERM(i > boxlo.x && i<boxhi.x, && j>boxlo.y && j<boxhi.y, && k>boxlo.z && k < boxhi.z))
                            {
                                divdw *= Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0);
                                Set::Vector gradpsi = Numeric::CellGradientOnNode(psi, i, j, k, 0, dx);
                                divdw += dw(i, j, k) * gradpsi;
                            }
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - divdw(d);
                        }
                    });
                }
            }

            a_rhs_mf[lev]->setMultiGhost(true);
            a_rhs_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
            a_ddw_mf[lev]->setMultiGhost(true);
            a_ddw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
        }
    }


public:
    Set::Scalar solve(const Set::Field<Set::Vector>& a_u_mf,
        const Set::Field<Set::Vector>& a_b_mf,
        Set::Field<T>& a_model_mf,
        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        if (verbose > 0)
            Util::Message(INFO, "Newton elastic operator: conservative_face_flux=",
                usesConservativeFaceFlux() ? 1 : 0);

        Set::Field<Set::Scalar> dsol_mf, rhs_mf, u0_mf, utrial_mf;
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;
        const int ddw_ncomp = usesConservativeFaceFlux() ? AMREX_SPACEDIM + 1 : 1;

        dsol_mf.resize(a_u_mf.finest_level + 1); dsol_mf.finest_level = a_u_mf.finest_level;
        dw_mf.resize(a_u_mf.finest_level + 1); dw_mf.finest_level = a_u_mf.finest_level;
        ddw_mf.resize(a_u_mf.finest_level + 1); ddw_mf.finest_level = a_u_mf.finest_level;
        rhs_mf.resize(a_u_mf.finest_level + 1); rhs_mf.finest_level = a_u_mf.finest_level;
        u0_mf.resize(a_u_mf.finest_level + 1); u0_mf.finest_level = a_u_mf.finest_level;
        utrial_mf.resize(a_u_mf.finest_level + 1); utrial_mf.finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            dsol_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf.NComp(),
                a_u_mf[lev]->nGrow());
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                ddw_ncomp,
                a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                a_b_mf.NComp(),
                a_b_mf[lev]->nGrow());
            u0_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf.NComp(),
                a_u_mf[lev]->nGrow());
            utrial_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf.NComp(),
                a_u_mf[lev]->nGrow());

            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            ddw_mf[lev]->setVal(Set::Matrix4<AMREX_SPACEDIM, T::sym>::Zero());

            a_b_mf.Copy(lev, *rhs_mf[lev], 0, 2);
            //amrex::MultiFab::Copy(*rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter + 1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);

            if (nriter == m_nriters) break;

            if (m_resync_coeffs && nriter > 0)
                m_elastic->SyncCoefficients();

            const Set::Scalar resnorm0 =
                m_line_search ? FieldNorm0(rhs_mf) : -1.0;

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                dsol_mf[lev]->setVal(0.0);
            Solver::Nonlocal::Linear::solve(dsol_mf, rhs_mf, a_tol_rel, a_tol_abs, checkpoint_file);

            Set::Scalar cornorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp, 0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;
                }
            }

            Set::Scalar alpha = 1.0;
            Set::Scalar resnorm = resnorm0;
            if (m_line_search)
            {
                for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                    a_u_mf.Copy(lev, *u0_mf[lev], 0, 2);

                auto restore_baseline = [&]()
                {
                    for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                        a_u_mf.CopyFrom(lev, *u0_mf[lev], 0, 2);
                };
                constexpr int max_backtrack = 8;
                bool accepted = false;
                for (int bt = 0; bt <= max_backtrack; bt++)
                {
                    for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                    {
                        amrex::MultiFab::Copy(*utrial_mf[lev], *u0_mf[lev],
                            0, 0, AMREX_SPACEDIM, 2);
                        amrex::MultiFab::Saxpy(*utrial_mf[lev], alpha,
                            *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
                        a_u_mf.CopyFrom(lev, *utrial_mf[lev], 0, 2);
                    }

                    try
                    {
                        prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf,
                            ddw_mf, a_model_mf);
                    }
                    catch (...)
                    {
                        restore_baseline();
                        prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf,
                            ddw_mf, a_model_mf);
                        throw;
                    }

                    resnorm = FieldNorm0(rhs_mf);

                    // Permit a small residual increase for roundoff-level
                    // changes in the composite AMR norm.
                    if (std::isfinite(resnorm0) && std::isfinite(resnorm) &&
                        (resnorm0 == 0.0 ? resnorm == 0.0 :
                            resnorm / resnorm0 <= 1.0001))
                    {
                        accepted = true;
                        break;
                    }

                    alpha *= 0.5;
                }
                if (!accepted)
                {
                    restore_baseline();
                    prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf,
                        ddw_mf, a_model_mf);
                    Util::Abort(INFO, "Newton line search failed after ",
                        max_backtrack, " backtracks: initial residual=", resnorm0,
                        ", final trial residual=", resnorm);
                }
            }
            else
            {
                for (int lev = 0; lev < dsol_mf.size(); ++lev)
                    a_u_mf.AddFrom(lev, *dsol_mf[lev], 0, 2);
                //amrex::MultiFab::Add(*a_u_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
            }

            const Set::Scalar full_update = cornorm;
            const Set::Scalar accepted_update = alpha * full_update;
            if (!std::isfinite(full_update) || !std::isfinite(accepted_update))
                Util::Abort(INFO, "Newton update norm is non-finite");
            if (verbose > 0)
            {
                Util::Message(INFO, "NR iteration ", nriter + 1,
                    ", alpha = ", alpha,
                    ", full max norm(ddisp) = ", full_update,
                    ", accepted max norm(ddisp) = ", accepted_update);
                if (m_line_search)
                    Util::Message(INFO, "NR residual: initial = ", resnorm0,
                        ", accepted = ", resnorm);
            }

            // Test the undamped Newton correction so backtracking cannot
            // manufacture convergence merely by choosing a very small alpha.
            if (m_nrtolerance > 0.0 && full_update < m_nrtolerance)
                return full_update;
            if (m_nrtolerance > 0.0 && nriter + 1 == m_nriters)
            {
                if (m_nr_abort_on_fail)
                    Util::Abort(INFO, "Newton solver failed to converge in ",
                        m_nriters, " iterations: final full update=", full_update,
                        ", final accepted update=", accepted_update,
                        ", tolerance=", m_nrtolerance);
                Util::Warning(INFO, "Newton solver failed to converge in ",
                    m_nriters, " iterations: final full update=", full_update,
                    ", final accepted update=", accepted_update,
                    ", tolerance=", m_nrtolerance);
                return full_update;
            }

        }

        return 0.0;
    }

    Set::Scalar solve(const Set::Field<Set::Scalar>& a_u_mf,
        const Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<T>& a_model_mf,
        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        Set::Field<Set::Scalar> dsol_mf, rhs_mf;
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;

        dsol_mf.resize(a_u_mf.finest_level + 1); dsol_mf.finest_level = a_u_mf.finest_level;
        dw_mf.resize(a_u_mf.finest_level + 1); dw_mf.finest_level = a_u_mf.finest_level;
        ddw_mf.resize(a_u_mf.finest_level + 1); ddw_mf.finest_level = a_u_mf.finest_level;
        rhs_mf.resize(a_u_mf.finest_level + 1); rhs_mf.finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            dsol_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf[lev]->nComp(),
                a_u_mf[lev]->nGrow());
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                a_b_mf[lev]->nComp(),
                a_b_mf[lev]->nGrow());

            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            ddw_mf[lev]->setVal(Set::Matrix4<AMREX_SPACEDIM, T::sym>::Zero());

            amrex::MultiFab::Copy(*rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter + 1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);


            if (nriter == m_nriters) break;

            Solver::Nonlocal::Linear::solve(dsol_mf, rhs_mf, a_tol_rel, a_tol_abs, checkpoint_file);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp, 0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;

                    Set::Scalar tmpsolnorm = a_u_mf[lev]->norm0(comp, 0);
                    if (tmpsolnorm > solnorm) solnorm = tmpsolnorm;
                }

            }
            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (verbose > 0) Util::Message(INFO, "NR iteration ", nriter + 1, ", relative norm(ddisp) = ", relnorm);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                amrex::MultiFab::Add(*a_u_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);

            if (relnorm < m_nrtolerance)
                return relnorm;

        }

        return 0.0;
    }
    Set::Scalar solve(const Set::Field<Set::Scalar>& a_u_mf,
        const Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<T>& a_model_mf)
    {
        return solve(a_u_mf, a_b_mf, a_model_mf, tol_rel, tol_abs);
    }

    void compResidual(Set::Field<Set::Scalar>& a_res_mf,
        Set::Field<Set::Scalar>& a_u_mf,
        Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<T>& a_model_mf)
    {
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;
        dw_mf.resize(a_u_mf.size());
        ddw_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }

        prepareForSolve(a_u_mf, a_b_mf, a_res_mf, dw_mf, ddw_mf, a_model_mf);
    }

    void compResidual(Set::Field<Set::Vector>& a_res_mf,
        Set::Field<Set::Vector>& a_u_mf,
        Set::Field<Set::Vector>& a_b_mf,
        Set::Field<T>& a_model_mf)
    {
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;
        Set::Field<Set::Scalar> res_mf(a_res_mf.size());
        dw_mf.resize(a_u_mf.size());
        ddw_mf.resize(a_u_mf.size());
        res_mf.resize(a_u_mf.size());
        const int ddw_ncomp = usesConservativeFaceFlux() ? AMREX_SPACEDIM + 1 : 1;
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                ddw_ncomp, a_b_mf[lev]->nGrow());
            res_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }

        prepareForSolve(a_u_mf, a_b_mf, res_mf, dw_mf, ddw_mf, a_model_mf);

        if (usesConservativeFaceFlux())
            for (int lev = a_u_mf.finest_level - 1; lev >= 0; --lev)
                m_elastic->Reflux(lev, *res_mf[lev], *res_mf[lev], *res_mf[lev],
                    *res_mf[lev + 1], *res_mf[lev + 1], *res_mf[lev + 1]);

        for (int lev = 0; lev < a_res_mf.size(); ++lev)
        {
            Util::RealFillBoundary(*res_mf[lev], m_elastic->Geom(lev));
            a_res_mf.CopyFrom(lev, *res_mf[lev], 0, 2);
        }
    }

    void compLinearSolverResidual(Set::Field<Set::Vector>& a_res_mf,
        Set::Field<Set::Vector>& a_u_mf,
        Set::Field<Set::Vector>& a_b_mf)
    {
        Util::Assert(INFO, TEST(a_res_mf.finest_level == a_u_mf.finest_level));
        Util::Assert(INFO, TEST(a_u_mf.finest_level == a_b_mf.finest_level));
        Set::Field<Set::Scalar> res_mf(a_res_mf.finest_level + 1), sol_mf(a_u_mf.finest_level + 1), rhs_mf(a_b_mf.finest_level + 1);
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            res_mf.Define(lev, a_b_mf[lev]->boxArray(), a_b_mf[lev]->DistributionMap(), AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            sol_mf.Define(lev, a_b_mf[lev]->boxArray(), a_b_mf[lev]->DistributionMap(), AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev, a_b_mf[lev]->boxArray(), a_b_mf[lev]->DistributionMap(), AMREX_SPACEDIM, a_b_mf[lev]->nGrow());

            a_u_mf.Copy(lev, *sol_mf[lev], 0, 2);
            a_b_mf.Copy(lev, *rhs_mf[lev], 0, 2);
        }

        mlmg->compResidual(amrex::GetVecOfPtrs(res_mf), amrex::GetVecOfPtrs(sol_mf), amrex::GetVecOfConstPtrs(rhs_mf));

        for (int lev = 0; lev <= a_res_mf.finest_level; ++lev)
        {
            a_res_mf.CopyFrom(lev, *res_mf[lev], 0, 2);
        }
    }


    void W(Set::Field<Set::Scalar>& a_w_mf,
        Set::Field<Set::Scalar>& a_u_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            BL_PROFILE("Solver::Nonlocal::Newton::DW()");

            const amrex::Real* DX = linop->Geom(lev).CellSize();
            amrex::Box domain(linop->Geom(lev).Domain());
            domain.convert(amrex::IntVect::TheNodeVector());

            for (MFIter mfi(*a_u_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                amrex::Array4<T> const& C = a_model_mf[lev]->array(mfi);
                amrex::Array4<amrex::Real> const& w = a_w_mf[lev]->array(mfi);
                amrex::Array4<const amrex::Real> const& u = a_u_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Matrix gradu;

                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    // Fill gradu
                    for (int p = 0; p < AMREX_SPACEDIM; p++)
                    {
                        AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX, sten));,
                            gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX, sten));,
                            gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX, sten)););
                    }

                    if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                        w(i, j, k) = C(i, j, k).W(gradu);
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        w(i, j, k) = C(i, j, k).W(0.5 * (gradu + gradu.transpose()));
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                        w(i, j, k) = C(i, j, k).W(gradu + Set::Matrix::Identity());
                });
            }
        }
    }

    void DW(Set::Field<Set::Scalar>& a_dw_mf,
        Set::Field<Set::Scalar>& a_u_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            BL_PROFILE("Solver::Nonlocal::Newton::DW()");

            const amrex::Real* DX = linop->Geom(lev).CellSize();
            amrex::Box domain(linop->Geom(lev).Domain());
            domain.convert(amrex::IntVect::TheNodeVector());

            for (MFIter mfi(*a_u_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                amrex::Array4<T> const& C = a_model_mf[lev]->array(mfi);
                amrex::Array4<amrex::Real> const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<const amrex::Real> const& u = a_u_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Matrix gradu;

                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    // Fill gradu
                    for (int p = 0; p < AMREX_SPACEDIM; p++)
                    {
                        AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX, sten));,
                            gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX, sten));,
                            gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX, sten)););
                    }

                    Set::Matrix sig = Set::Matrix::Zero();

                    if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                        sig = C(i, j, k).DW(gradu);
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        sig = C(i, j, k).DW(0.5 * (gradu + gradu.transpose()));
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                        sig = C(i, j, k).DW(gradu + Set::Matrix::Identity());

                    // = C(i,j,k)(gradu,m_homogeneous);

                    AMREX_D_PICK(dw(i, j, k, 0) = sig(0, 0);
                    ,
                        dw(i, j, k, 0) = sig(0, 0); dw(i, j, k, 1) = sig(0, 1);
                    dw(i, j, k, 2) = sig(1, 0); dw(i, j, k, 3) = sig(1, 1);
                    ,
                        dw(i, j, k, 0) = sig(0, 0); dw(i, j, k, 1) = sig(0, 1); dw(i, j, k, 2) = sig(0, 2);
                    dw(i, j, k, 3) = sig(1, 0); dw(i, j, k, 4) = sig(1, 1); dw(i, j, k, 5) = sig(1, 2);
                    dw(i, j, k, 6) = sig(2, 0); dw(i, j, k, 7) = sig(2, 1); dw(i, j, k, 8) = sig(2, 2););

                });
            }
        }
    }


public:
    int m_nriters = 1;
    Set::Scalar m_nrtolerance = 0.0;
    bool m_line_search = false;
    bool m_resync_coeffs = false;
    bool m_nr_abort_on_fail = false;
    Operator::Elastic<T::sym>* m_elastic;
    //BC::Operator::Elastic::Elastic *m_bc;

    Set::Field<Set::Scalar>* m_psi = nullptr;
    bool m_conservative_face_flux = false;

private:
    bool m_reported_psi_face_flux_conflict = false;

    Set::Scalar FieldNorm0(const Set::Field<Set::Scalar>& a_mf) const
    {
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> composite(a_mf.size());
        for (int lev = 0; lev < a_mf.size(); ++lev)
        {
            composite[lev] = std::make_unique<amrex::MultiFab>(
                a_mf[lev]->boxArray(), a_mf[lev]->DistributionMap(),
                a_mf[lev]->nComp(), a_mf[lev]->nGrowVect());
            amrex::MultiFab::Copy(*composite[lev], *a_mf[lev], 0, 0,
                a_mf[lev]->nComp(), a_mf[lev]->nGrowVect());
        }

        for (int lev = static_cast<int>(composite.size()) - 2; lev >= 0; --lev)
            m_elastic->Reflux(lev, *composite[lev], *composite[lev],
                *composite[lev], *composite[lev + 1], *composite[lev + 1],
                *composite[lev + 1]);

        Set::Scalar norm = 0.0;
        for (int lev = 0; lev < composite.size(); ++lev)
        {
            const Set::Scalar tmp =
                m_elastic->normInf(lev, *composite[lev], false);
            if (tmp > norm) norm = tmp;
        }
        return norm;
    }

public:
    // These paramters control a standard Newton-Raphson solve.
    // 
    // **Note**: 
    // This class inherits all of the linear solve paramters
    // from its parent class, :ref:`Solver::Nonlocal::Linear`
    static void Parse(Newton<T>& value, amrex::ParmParse& pp)
    {
        Linear::Parse(value, pp);

        // Number of newton-raphson iterations.
        pp_query("nriters", value.m_nriters);

        // Tolerance for the undamped absolute correction norm in vector
        // mechanics solves, and the relative correction norm in scalar solves.
        pp_query("nrtolerance", value.m_nrtolerance); 

        // Backtrack nonlinear updates until the residual is within a small
        // relative allowance of its initial value.
        pp_query("line_search", value.m_line_search);

        // Rebuild coarse MG coefficients and the smoother diagonal after each
        // nonlinear relinearization. MLMG otherwise prepares them only once.
        pp_query("resync_coeffs", value.m_resync_coeffs);

        // Abort instead of warning if Newton iterations do not converge.
        pp_query("nr_abort_on_fail", value.m_nr_abort_on_fail);

        // Discretize the elastic residual with coefficient-aware face fluxes.
        pp_query("conservative_face_flux", value.m_conservative_face_flux);
    }

};
} // namespace Nonlocal
} // namespace Solver


#endif
