#include "Elastic.H"
#include "AMReX_Loop.H"
#include "Set/Set.H"

#include "Numeric/Stencil.H"

#ifdef ALAMO_GPU
#define ALAMO_ELASTIC_OP_FOR amrex::ParallelFor
#define ALAMO_ELASTIC_OP_CAPTURE [=]
#define ALAMO_ELASTIC_OP_DEVICE AMREX_GPU_DEVICE
#define ALAMO_ELASTIC_OP_BC_EVAL(bc, bc_type, u, gradu, sigma, i, j, k, bx) \
    ::BC::Operator::Elastic::Elastic::eval(bc_type, u, gradu, sigma, i, j, k, bx)
#else
#define ALAMO_ELASTIC_OP_FOR amrex::LoopConcurrentOnCpu
#define ALAMO_ELASTIC_OP_CAPTURE [=]
#define ALAMO_ELASTIC_OP_DEVICE
#define ALAMO_ELASTIC_OP_BC_EVAL(bc, bc_type, u, gradu, sigma, i, j, k, bx) \
    (*(bc))(u, gradu, sigma, i, j, k, bx)
#endif

namespace Operator
{
template<int SYM>
Elastic<SYM>::Elastic(const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info)
{
    BL_PROFILE("Operator::Elastic::Elastic()");

    define(a_geom, a_grids, a_dmap, a_info);
}

template<int SYM>
Elastic<SYM>::~Elastic()
{}

template<int SYM>
void
Elastic<SYM>::define(const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("Operator::Elastic::define()");

    Operator::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    int model_nghost = 2;

    m_ddw_mf.resize(m_num_amr_levels);
    m_psi_mf.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_ddw_mf[amrlev].resize(m_num_mg_levels[amrlev]);
        m_psi_mf[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_ddw_mf[amrlev][mglev].reset(new MultiTab(amrex::convert(m_grids[amrlev][mglev],
                amrex::IntVect::TheNodeVector()),
                m_dmap[amrlev][mglev], AMREX_SPACEDIM + 1, model_nghost));
            m_psi_mf[amrlev][mglev].reset(new MultiFab(m_grids[amrlev][mglev],
                m_dmap[amrlev][mglev], 1, model_nghost));

            if (!m_psi_set) m_psi_mf[amrlev][mglev]->setVal(1.0);
        }
    }
}

template <int SYM>
void
Elastic<SYM>::SetModel(MATRIX4& a_model)
{
    for (int amrlev = 0; amrlev < m_num_amr_levels; amrlev++)
    {
        amrex::Box domain(m_geom[amrlev][0].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::DeviceErrorFlag setmodel_error;
        int* setmodel_error_flag = setmodel_error.dataPtr();
#endif
#endif
        for (MFIter mfi(*m_ddw_mf[amrlev][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.grownnodaltilebox();

            amrex::Array4<MATRIX4> const& ddw = (*(m_ddw_mf[amrlev][0])).array(mfi);

            ALAMO_ELASTIC_OP_FOR(bx, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k) {
                for (int n = 0; n < AMREX_SPACEDIM + 1; ++n)
                    ddw(i, j, k, n) = a_model;

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
                if (ddw(i, j, k, 0).contains_nan()) Util::SetDeviceError(setmodel_error_flag);
#else
                if (ddw(i, j, k, 0).contains_nan()) Util::Abort(INFO, "model is nan at (", i, ",", j, ",", k, "), amrlev=", amrlev);
#endif
#endif
            });
        }
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(setmodel_error, INFO, "Operator::Elastic::SetModel() detected a NaN in the device coefficient field");
#endif
#endif
        (*(m_ddw_mf[amrlev][0])).setMultiGhost(true);
        (*(m_ddw_mf[amrlev][0])).FillBoundary( Geom(amrlev,0).periodicity());
    }
    m_model_set = true;
}

template <int SYM>
void
Elastic<SYM>::SetModel(int amrlev, const amrex::FabArray<amrex::BaseFab<MATRIX4> >& a_model)
{
    BL_PROFILE("Operator::Elastic::SetModel()");

    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());

    if (a_model.boxArray() != m_ddw_mf[amrlev][0]->boxArray()) Util::Abort(INFO, "Inconsistent box arrays\n", "a_model.boxArray()=\n", a_model.boxArray(), "\n but the current box array is \n", m_ddw_mf[amrlev][0]->boxArray());
    if (a_model.DistributionMap() != m_ddw_mf[amrlev][0]->DistributionMap()) Util::Abort(INFO, "Inconsistent distribution maps");
    if (a_model.nComp() != 1 &&
        a_model.nComp() != m_ddw_mf[amrlev][0]->nComp())
        Util::Abort(INFO, "Inconsistent # of coefficient components - should be 1 or ",
            m_ddw_mf[amrlev][0]->nComp());
    if (a_model.nGrow() != m_ddw_mf[amrlev][0]->nGrow()) Util::Abort(INFO, "Inconsistent # of ghost nodes, should be ", m_ddw_mf[amrlev][0]->nGrow());

    const bool nodal_only = a_model.nComp() == 1;

    for (MFIter mfi(a_model, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.grownnodaltilebox();

        amrex::Array4<MATRIX4> const& C = (*(m_ddw_mf[amrlev][0])).array(mfi);
        amrex::Array4<const MATRIX4> const& a_C = a_model.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            C(i, j, k, 0) = a_C(i, j, k, 0);
            for (int n = 1; n < AMREX_SPACEDIM + 1; ++n)
                C(i, j, k, n) = a_C(i, j, k, nodal_only ? 0 : n);
        });
    }
    m_ddw_mf[amrlev][0]->setMultiGhost(true);
    m_ddw_mf[amrlev][0]->FillBoundaryAndSync(Geom(amrlev,0).periodicity());
    m_model_set = true;
}

template <int SYM>
void
Elastic<SYM>::SetPsi(int amrlev, const amrex::MultiFab& a_psi_mf)
{
    BL_PROFILE("Operator::Elastic::SetPsi()");
    amrex::Box domain(m_geom[amrlev][0].Domain());

    for (MFIter mfi(a_psi_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox() & domain;

        amrex::Array4<Set::Scalar> const& m_psi = (*(m_psi_mf[amrlev][0])).array(mfi);
        amrex::Array4<const Set::Scalar> const& a_psi = a_psi_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            m_psi(i, j, k) = a_psi(i, j, k);
        });
    }
    m_psi_set = true;
}

template<int SYM>
void
Elastic<SYM>::Fapply(int amrlev, int mglev, MultiFab& a_f, const MultiFab& a_u) const
{
    BL_PROFILE("Operator::Elastic::Fapply()");

    amrex::Box domain(m_geom[amrlev][mglev].growPeriodicDomain(1));
    domain.convert(amrex::IntVect::TheNodeVector());

    amrex::Box stencilbox(m_geom[amrlev][mglev].growPeriodicDomain(2));
    stencilbox.convert(amrex::IntVect::TheNodeVector());

    Set::Vector DX(m_geom[amrlev][mglev].CellSize()); // device-safe cell size (copied into device closure)
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
    Util::DeviceErrorFlag fapply_error;
    int* fapply_error_flag = fapply_error.dataPtr();
#endif
#endif

    for (MFIter mfi(a_f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox().grow(1) & domain;
        amrex::Box tilebox = mfi.grownnodaltilebox() & bx;

        amrex::Array4<MATRIX4> const& DDW = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
        amrex::Array4<const amrex::Real> const& U = a_u.array(mfi);
        amrex::Array4<amrex::Real> const& F = a_f.array(mfi);
        amrex::Array4<Set::Scalar> const& psi = m_psi_mf[amrlev][mglev]->array(mfi);

        const Dim3 lo = amrex::lbound(stencilbox), hi = amrex::ubound(stencilbox);
        auto m_psi_set = this->m_psi_set;
        auto m_psi_small = this->m_psi_small;
        auto m_conservative_face_flux = this->m_conservative_face_flux;
        auto m_uniform = this->m_uniform;
#ifdef ALAMO_GPU
        auto m_bc_type = this->m_bc->GetBcTypeArray();
#else
        auto m_bc = this->m_bc;
#endif

        ALAMO_ELASTIC_OP_FOR(tilebox, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k) {

            Set::Vector f = Set::Vector::Zero();

            Set::Vector u;
            for (int p = 0; p < AMREX_SPACEDIM; p++) u(p) = U(i, j, k, p);


            bool AMREX_D_DECL(xmin = (i == lo.x), ymin = (j == lo.y), zmin = (k == lo.z)),
                AMREX_D_DECL(xmax = (i == hi.x), ymax = (j == hi.y), zmax = (k == hi.z));

            std::array<Numeric::StencilType, AMREX_SPACEDIM>
                sten = Numeric::GetStencil(i, j, k, stencilbox);

            Set::Matrix gradu; // gradu(i,j) = u_{i,j)

            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(U, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(U, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(U, i, j, k, p, DX.data(), sten)););
            }

            Set::Scalar psi_avg = 1.0;
            if (m_psi_set) psi_avg = (1.0 - m_psi_small) * Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0) + m_psi_small;

            // ddw is reused below (sig, C(gradgradu), grad(psi) correction) - loaded once.
            MATRIX4 const ddw = DDW(i, j, k);

            amrex::IntVect m(AMREX_D_DECL(i, j, k));
            if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
            {
                // Only boundary rows consume the stress tensor.
                Set::Matrix sig = (ddw * gradu) * psi_avg;
                f = ALAMO_ELASTIC_OP_BC_EVAL(m_bc, m_bc_type, u, gradu, sig, i, j, k, stencilbox);
            }
            else
            {
                if (m_conservative_face_flux)
                {
                    // Centered tangential endpoint averages avoid the checkerboard nullspace.
                    for (int face = 0; face < AMREX_SPACEDIM; ++face)
                    {
                        const int im = i - (face == 0);
                        const int jm = j - (face == 1);
                        const int km = k - (face == 2);
                        const Set::Matrix grad_hi = Numeric::FaceGradient(
                            U, i, j, k, face, DX.data());
                        const Set::Matrix grad_lo = Numeric::FaceGradient(
                            U, im, jm, km, face, DX.data());
                        const Set::Matrix flux_hi =
                            DDW(i, j, k, face + 1) * grad_hi;
                        const Set::Matrix flux_lo =
                            DDW(im, jm, km, face + 1) * grad_lo;
                        f += (flux_hi.col(face) - flux_lo.col(face)) / DX[face];
                    }
                }
                else
                {

                /// \todo replace with Numeric::Hessian(U,i,j,k,DX,sten)
                Set::Matrix3 gradgradu; // gradgradu[k](l,j) = u_{k,lj}

                for (int p = 0; p < AMREX_SPACEDIM; p++)
                {
                    AMREX_D_TERM(
                        gradgradu(p, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(U, i, j, k, p, DX.data()));,
                        gradgradu(p, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(U, i, j, k, p, DX.data()));,
                        gradgradu(p, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(U, i, j, k, p, DX.data())););

                    AMREX_D_TERM(
                        ,
                        gradgradu(p, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(U, i, j, k, p, DX.data()));
                        gradgradu(p, 1, 0) = gradgradu(p, 0, 1);
                        ,
                        gradgradu(p, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(U, i, j, k, p, DX.data()));
                        gradgradu(p, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(U, i, j, k, p, DX.data()));
                        gradgradu(p, 2, 0) = gradgradu(p, 0, 2);
                        gradgradu(p, 2, 1) = gradgradu(p, 1, 2););
                }

                // f = C(grad grad u) + grad(C)*grad(u)
                f = (ddw * gradgradu) * psi_avg;

                if (!m_uniform)
                {
                    // Accumulate grad(C):grad(u) one spatial direction at a time so
                    // that only a single Matrix4 derivative temp is live at any
                    // moment. Naming three temps (Cgrad1/2/3 = 135 live doubles in
                    // 3D) whose scope spanned the whole expression was the dominant
                    // Fapply register-spill source (255 regs/thread -> ~12.5%
                    // occupancy, PHASE_A_FINDINGS.md sec.4).
                    //
                    // MulCol(a,b,c) is bit-identical to (a*b).col(c), computed
                    // directly, and the summation order is unchanged, so f stays
                    // bit-identical to both prior forms.
                    Set::Vector graddc = Set::Vector::Zero();
                    graddc += Set::MulCol(Numeric::Stencil<MATRIX4, 1, 0, 0>::D(DDW, i, j, k, 0, DX.data(), sten), gradu, 0);
#if AMREX_SPACEDIM > 1
                    graddc += Set::MulCol(Numeric::Stencil<MATRIX4, 0, 1, 0>::D(DDW, i, j, k, 0, DX.data(), sten), gradu, 1);
#endif
#if AMREX_SPACEDIM > 2
                    graddc += Set::MulCol(Numeric::Stencil<MATRIX4, 0, 0, 1>::D(DDW, i, j, k, 0, DX.data(), sten), gradu, 2);
#endif
                    f += graddc * psi_avg;
                }
                if (m_psi_set)
                {
                    Set::Vector gradpsi = Numeric::CellGradientOnNode(psi, i, j, k, 0, DX.data());
                    gradpsi *= (1.0 - m_psi_small);
                    f += (ddw * gradu) * gradpsi;
                }
                }

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
                if (std::isnan(f(0)) || std::isnan(f(1))
#if AMREX_SPACEDIM == 3
                    || std::isnan(f(2))
#endif
                )
                {
                    Util::SetDeviceError(fapply_error_flag);
                }
#else
                if (std::isnan(f(0)) || std::isnan(f(1)))
                {
                    Util::Message(INFO,"  =================  ");
                    Util::Message(INFO,"amrlev=",amrlev);
                    Util::Message(INFO,"mglev=",mglev);
                    Util::Message(INFO,"i=",i," j=",j);
                    Util::Message(INFO,"f:            ",f.transpose());
                    Util::Message(INFO,"U(i,j,k):     ",U(i,j,k,0)," ",U(i,j,k,1));
                    Util::Message(INFO,"U(i-1,j,k):   ",U(i-1,j,k,0)," ",U(i-1,j,k,1));
                    Util::Message(INFO,"U(i+1,j,k):   ",U(i+1,j,k,0)," ",U(i-1,j,k,1));
                    Util::Message(INFO,"U(i,j-1,k):     ",U(i,j-1,k,0)," ",U(i,j-1,k,1));
                    Util::Message(INFO,"U(i,j+1,k):     ",U(i,j+1,k,0)," ",U(i,j+1,k,1));
                    Util::Message(INFO,"gradu:        ",gradu);
                    Util::Message(INFO,"gradgradu[0]: ",gradgradu[0]);
                    Util::Message(INFO,"gradgradu[1]: ",gradgradu[1]);
                    Util::Message(INFO,"DDW (i  ,j  ): ",DDW(i,j,k));
                    Util::Message(INFO,"DDW (i-1,j  ): ",DDW(i-1,j,k));
                    Util::Message(INFO,"DDW (i+1,j  ): ",DDW(i+1,j,k));
                    Util::Message(INFO,"DDW (i  ,j+1): ",DDW(i,j+1,k));
                    Util::Message(INFO,"DDW (i  ,j-1): ",DDW(i,j-1,k));
                    Util::Message(INFO,"psi_av: ",psi_avg);
                    Util::Message(INFO,"psi_set: ",m_psi_set);
                    Util::Message(INFO,"  =================  ");
                    Util::Abort(INFO);
                }
#endif
#endif
            }
            AMREX_D_TERM(F(i, j, k, 0) = f[0];, F(i, j, k, 1) = f[1];, F(i, j, k, 2) = f[2];);
        });
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(fapply_error, INFO, "Operator::Elastic::Fapply() detected an invalid value on device");
#endif
#endif
    }

}



template<int SYM>
void
Elastic<SYM>::Diagonal(int amrlev, int mglev, MultiFab& a_diag)
{
    BL_PROFILE("Operator::Elastic::Diagonal()");

    const amrex::IntVect diagonal_nghost = a_diag.nGrowVect();
    amrex::Box domain(m_geom[amrlev][mglev].growPeriodicDomain(
        diagonal_nghost.max()));
    domain.convert(amrex::IntVect::TheNodeVector());

    amrex::Box stencilbox(m_geom[amrlev][mglev].growPeriodicDomain(2));
    stencilbox.convert(amrex::IntVect::TheNodeVector());

    Set::Vector DX(m_geom[amrlev][mglev].CellSize()); // device-safe cell size (copied into device closure)
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
    Util::DeviceErrorFlag diagonal_error;
    int* diagonal_error_flag = diagonal_error.dataPtr();
#endif
#endif

    for (MFIter mfi(a_diag, false); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox().grow(diagonal_nghost) & domain;
        amrex::Box tilebox = mfi.grownnodaltilebox() & bx;

        amrex::Array4<MATRIX4> const& DDW = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
        amrex::Array4<Set::Scalar> const& diag = a_diag.array(mfi);
        amrex::Array4<Set::Scalar> const& psi = m_psi_mf[amrlev][mglev]->array(mfi);

        const Dim3 lo = amrex::lbound(stencilbox), hi = amrex::ubound(stencilbox);
        auto m_psi_set = this->m_psi_set;
        auto m_psi_small = this->m_psi_small;
        auto m_conservative_face_flux = this->m_conservative_face_flux;
#ifdef ALAMO_GPU
        auto m_bc_type = this->m_bc->GetBcTypeArray();
#else
        auto m_bc = this->m_bc;
#endif

        amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            std::array<Numeric::StencilType, AMREX_SPACEDIM>
                sten = Numeric::GetStencil(i, j, k, stencilbox);

            // gradu(i,j) = u_{i,j)
            std::array<Set::Matrix,AMREX_SPACEDIM> gradu = Numeric::Gradient_Diagonal<Set::Matrix>(DX.data(), sten);

            // gradgradu[k](l,j) = u_{k,lj}
            std::array<Set::Matrix3,AMREX_SPACEDIM>  gradgradu = Numeric::Gradient_Diagonal<Set::Matrix3>(DX.data());


            Set::Vector f = Set::Vector::Zero();

            bool
                AMREX_D_DECL(xmin = (i == lo.x), ymin = (j == lo.y), zmin = (k == lo.z)),
                AMREX_D_DECL(xmax = (i == hi.x), ymax = (j == hi.y), zmax = (k == hi.z));

            Set::Scalar psi_avg = 1.0;
            if (m_psi_set) psi_avg = (1.0 - m_psi_small) * Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0) + m_psi_small;

            MATRIX4 const ddw = DDW(i, j, k);

            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {

                diag(i, j, k, p) = 0.0;


                amrex::IntVect m(AMREX_D_DECL(i, j, k));
                if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
                {
                    Set::Matrix sig = ddw * gradu[p] * psi_avg;
                    Set::Vector u = Set::Vector::Zero();
                    u(p) = 1.0;
                    f = ALAMO_ELASTIC_OP_BC_EVAL(m_bc, m_bc_type, u, gradu[p], sig, i, j, k, stencilbox);
                    diag(i, j, k, p) = f(p);
                }
                else if (m_conservative_face_flux)
                {
                    // Only the two normal face differences contribute to this diagonal entry.
                    for (int face = 0; face < AMREX_SPACEDIM; ++face)
                    {
                        const int im = i - (face == 0);
                        const int jm = j - (face == 1);
                        const int km = k - (face == 2);
                        diag(i, j, k, p) -=
                            (DDW(i, j, k, face + 1)(p, face, p, face)
                             + DDW(im, jm, km, face + 1)(p, face, p, face))
                            / (DX[face] * DX[face]);
                    }
                }
                else
                {
                    Set::Vector f = (ddw * gradgradu[p]) * psi_avg;
                    diag(i, j, k, p) += f(p);
                }

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
                if (std::isnan(diag(i, j, k, p)) || std::isinf(diag(i, j, k, p)) || diag(i, j, k, p) == 0)
                {
                    Util::SetDeviceError(diagonal_error_flag);
                }
#else
                if (std::isnan(diag(i, j, k, p))) Util::Abort(INFO, "diagonal is nan at (", i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
                if (std::isinf(diag(i, j, k, p))) Util::Abort(INFO, "diagonal is inf at (", i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
                if (diag(i, j, k, p) == 0) Util::Abort(INFO, "diagonal is zero at (", i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
#endif
#endif

            }
        });
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(diagonal_error, INFO, "Operator::Elastic::Diagonal() detected an invalid value on device");
#endif
#endif
    }

    a_diag.FillBoundaryAndSync(Geom(amrlev,mglev).periodicity());
    nodalSync(amrlev,mglev,a_diag);
}


template<int SYM>
void
Elastic<SYM>::Error0x(int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const
{
    BL_PROFILE("Operator::Elastic::Error0x()");
    Util::Message(INFO);

    int ncomp = x.nComp();//getNComp();
    int nghost = x.nGrow();

    if (!m_diagonal_computed)
        Util::Abort(INFO, "Operator::Diagonal() must be called before using normalize");

    amrex::MultiFab D0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);
    amrex::MultiFab AD0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);

    amrex::MultiFab::Copy(D0x, x, 0, 0, ncomp, nghost); // D0x = x
    amrex::MultiFab::Divide(D0x, *m_diag[amrlev][mglev], 0, 0, ncomp, 0); // D0x = x/diag
    amrex::MultiFab::Copy(AD0x, D0x, 0, 0, ncomp, nghost); // AD0x = D0x

    Fapply(amrlev, mglev, AD0x, D0x);    // AD0x = A * D0 * x

    amrex::MultiFab::Copy(R0x, x, 0, 0, ncomp, nghost); // R0x = x
    amrex::MultiFab::Subtract(R0x, AD0x, 0, 0, ncomp, nghost); // R0x = x - AD0x
}


template<int SYM>
void
Elastic<SYM>::FFlux(int /*amrlev*/, const MFIter& /*mfi*/,
    const std::array<FArrayBox*, AMREX_SPACEDIM>& sigmafab,
    const FArrayBox& /*ufab*/, const int /*face_only*/) const
{
    BL_PROFILE("Operator::Elastic::FFlux()");
    Util::Message(INFO);
    amrex::BaseFab<amrex::Real> AMREX_D_DECL(&fxfab = *sigmafab[0],
        &fyfab = *sigmafab[1],
        &fzfab = *sigmafab[2]);
    AMREX_D_TERM(fxfab.setVal<amrex::RunOn::Device>(0.0);,
        fyfab.setVal<amrex::RunOn::Device>(0.0);,
        fzfab.setVal<amrex::RunOn::Device>(0.0););

}

template<int SYM>
void
Elastic<SYM>::Strain(int amrlev,
    amrex::MultiFab& a_eps,
    const amrex::MultiFab& a_u,
    bool voigt) const
{
    BL_PROFILE("Operator::Elastic::Strain()");

    Set::Vector DX(m_geom[amrlev][0].CellSize()); // device-safe cell size (copied into device closure)
    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());


    for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& epsilon = a_eps.array(mfi);
        amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Matrix gradu;

            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten
                = Numeric::GetStencil(i, j, k, domain);

            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
            }

            Set::Matrix eps = 0.5 * (gradu + gradu.transpose());

            if (voigt)
            {
                AMREX_D_PICK(epsilon(i, j, k, 0) = eps(0, 0);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(1, 1); epsilon(i, j, k, 2) = eps(0, 1);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(1, 1); epsilon(i, j, k, 2) = eps(2, 2);
                epsilon(i, j, k, 3) = eps(1, 2); epsilon(i, j, k, 4) = eps(2, 0); epsilon(i, j, k, 5) = eps(0, 1););
            }
            else
            {
                AMREX_D_PICK(epsilon(i, j, k, 0) = eps(0, 0);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(0, 1);
                epsilon(i, j, k, 2) = eps(1, 0); epsilon(i, j, k, 3) = eps(1, 1);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(0, 1); epsilon(i, j, k, 2) = eps(0, 2);
                epsilon(i, j, k, 3) = eps(1, 0); epsilon(i, j, k, 4) = eps(1, 1); epsilon(i, j, k, 5) = eps(1, 2);
                epsilon(i, j, k, 6) = eps(2, 0); epsilon(i, j, k, 7) = eps(2, 1); epsilon(i, j, k, 8) = eps(2, 2););
            }
        });
    }
}


template<int SYM>
void
Elastic<SYM>::Stress(int amrlev,
    amrex::MultiFab& a_sigma,
    const amrex::MultiFab& a_u,
    bool voigt, bool a_homogeneous)
{
    BL_PROFILE("Operator::Elastic::Stress()");
    SetHomogeneous(a_homogeneous);

    Set::Vector DX(m_geom[amrlev][0].CellSize()); // device-safe cell size (copied into device closure)
    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());

    for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& DDW = (*(m_ddw_mf[amrlev][0])).array(mfi);
        amrex::Array4<amrex::Real> const& sigma = a_sigma.array(mfi);
        amrex::Array4<Set::Scalar> const& psi = m_psi_mf[amrlev][0]->array(mfi);
        amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
        auto m_psi_set = this->m_psi_set;
        auto m_psi_small = this->m_psi_small;
        ALAMO_ELASTIC_OP_FOR(bx, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k)
        {
            Set::Matrix gradu;

            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten
                = Numeric::GetStencil(i, j, k, domain);

            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
            }

            Set::Scalar psi_avg = 1.0;
            if (m_psi_set) psi_avg = (1.0 - m_psi_small) * Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0) + m_psi_small;
            Set::Matrix sig = (DDW(i, j, k) * gradu) * psi_avg;

            if (voigt)
            {
                AMREX_D_PICK(sigma(i, j, k, 0) = sig(0, 0);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(1, 1); sigma(i, j, k, 2) = sig(0, 1);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(1, 1); sigma(i, j, k, 2) = sig(2, 2);
                sigma(i, j, k, 3) = sig(1, 2); sigma(i, j, k, 4) = sig(2, 0); sigma(i, j, k, 5) = sig(0, 1););
            }
            else
            {
                AMREX_D_PICK(sigma(i, j, k, 0) = sig(0, 0);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(0, 1);
                sigma(i, j, k, 2) = sig(1, 0); sigma(i, j, k, 3) = sig(1, 1);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(0, 1); sigma(i, j, k, 2) = sig(0, 2);
                sigma(i, j, k, 3) = sig(1, 0); sigma(i, j, k, 4) = sig(1, 1); sigma(i, j, k, 5) = sig(1, 2);
                sigma(i, j, k, 6) = sig(2, 0); sigma(i, j, k, 7) = sig(2, 1); sigma(i, j, k, 8) = sig(2, 2););
            }
        });
    }
}


template<int SYM>
void
Elastic<SYM>::Energy(int amrlev,
    amrex::MultiFab& a_energy,
    const amrex::MultiFab& a_u, bool a_homogeneous)
{
    BL_PROFILE("Operator::Elastic::Energy()");
    SetHomogeneous(a_homogeneous);

    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());

    Set::Vector DX(m_geom[amrlev][0].CellSize()); // device-safe cell size (copied into device closure)

    for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& DDW = (*(m_ddw_mf[amrlev][0])).array(mfi);
        amrex::Array4<amrex::Real> const& energy = a_energy.array(mfi);
        amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Matrix gradu;

            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten
                = Numeric::GetStencil(i, j, k, domain);

            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
            }

            Set::Matrix eps = .5 * (gradu + gradu.transpose());
            Set::Matrix sig = DDW(i, j, k) * gradu;

            for (int m = 0; m < AMREX_SPACEDIM; m++)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++)
                {
                    energy(i, j, k) += .5 * sig(m, n) * eps(m, n);
                }
            }
        });
    }
}

template<int SYM>
void
Elastic<SYM>::averageDownCoeffs()
{
    BL_PROFILE("Elastic::averageDownCoeffs()");

    if (m_average_down_coeffs)
        for (int amrlev = m_num_amr_levels - 1; amrlev > 0; --amrlev)
            averageDownCoeffsDifferentAmrLevels(amrlev);

    averageDownCoeffsSameAmrLevel(0);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            if (m_ddw_mf[amrlev][mglev]) {
                FillBoundaryCoeff(*m_ddw_mf[amrlev][mglev], Geom(amrlev,mglev).periodicity());
                FillBoundaryCoeff(*m_psi_mf[amrlev][mglev], Geom(amrlev,mglev).periodicity());
            }
        }
    }

}

template<int SYM>
void
Elastic<SYM>::averageDownCoeffsDifferentAmrLevels(int fine_amrlev)
{
    BL_PROFILE("Operator::Elastic::averageDownCoeffsDifferentAmrLevels()");
    Util::Assert(INFO, TEST(fine_amrlev > 0));

    const int crse_amrlev = fine_amrlev - 1;
    MultiTab& crse_ddw = *m_ddw_mf[crse_amrlev][0];
    MultiTab& fine_ddw = *m_ddw_mf[fine_amrlev][0];
    const int ncomp = crse_ddw.nComp();

    amrex::Box cdomain(m_geom[crse_amrlev][0].Domain());
    cdomain.convert(amrex::IntVect::TheNodeVector());

    const Geometry& cgeom = m_geom[crse_amrlev][0];

    const BoxArray& fba = fine_ddw.boxArray();
    const DistributionMapping& fdm = fine_ddw.DistributionMap();

    MultiTab fine_ddw_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
    fine_ddw_for_coarse.ParallelCopy(crse_ddw, 0, 0, ncomp, 0, 0, cgeom.periodicity());

    const int coarse_fine_node = 1;
    const int fine_fine_node = 2;

    amrex::iMultiFab nodemask(amrex::coarsen(fba, 2), fdm, 1, 2);
    nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev], 0, 0, 1, 0, 0, cgeom.periodicity());

    amrex::iMultiFab cellmask(amrex::convert(amrex::coarsen(fba, 2), amrex::IntVect::TheCellVector()), fdm, 1, 2);
    cellmask.ParallelCopy(*m_cc_fine_mask[crse_amrlev], 0, 0, 1, 1, 1, cgeom.periodicity());

    for (MFIter mfi(fine_ddw_for_coarse, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        amrex::Array4<const int> const& nmask = nodemask.array(mfi);

        amrex::Array4<MATRIX4> const& cdata = fine_ddw_for_coarse.array(mfi);
        amrex::Array4<const MATRIX4> const& fdata = fine_ddw.array(mfi);

        const Dim3 lo = amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

        for (int n = 0; n < ncomp; n++)
        {
            // (I,J,K) = coarse, (i,j,k) = fine
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int I, int J, int K) {
                int i = I * 2, j = J * 2, k = K * 2;

                if (nmask(I, J, K) == fine_fine_node || nmask(I, J, K) == coarse_fine_node)
                {
                    if (n > 0)
                    {
                        const int face = n - 1;
                        cdata(I, J, K, n) = 0.5 * (
                            fdata(i, j, k, n)
                            + fdata(i + (face == 0), j + (face == 1),
                                k + (face == 2), n));
                        return;
                    }
                    if ((I == lo.x || I == hi.x) &&
                        (J == lo.y || J == hi.y) &&
                        (K == lo.z || K == hi.z)) // Corner
                        cdata(I, J, K, n) = fdata(i, j, k, n);
                    else if ((J == lo.y || J == hi.y) &&
                        (K == lo.z || K == hi.z)) // X edge
                        cdata(I, J, K, n) = fdata(i - 1, j, k, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i + 1, j, k, n) * 0.25;
                    else if ((K == lo.z || K == hi.z) &&
                        (I == lo.x || I == hi.x)) // Y edge
                        cdata(I, J, K, n) = fdata(i, j - 1, k, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i, j + 1, k, n) * 0.25;
                    else if ((I == lo.x || I == hi.x) &&
                        (J == lo.y || J == hi.y)) // Z edge
                        cdata(I, J, K, n) = fdata(i, j, k - 1, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i, j, k + 1, n) * 0.25;
                    else if (I == lo.x || I == hi.x) // X face
                        cdata(I, J, K, n) =
                        (fdata(i, j - 1, k - 1, n) + fdata(i, j, k - 1, n) * 2.0 + fdata(i, j + 1, k - 1, n)
                            + fdata(i, j - 1, k, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i, j + 1, k, n) * 2.0
                            + fdata(i, j - 1, k + 1, n) + fdata(i, j, k + 1, n) * 2.0 + fdata(i, j + 1, k + 1, n)) / 16.0;
                    else if (J == lo.y || J == hi.y) // Y face
                        cdata(I, J, K, n) =
                        (fdata(i - 1, j, k - 1, n) + fdata(i - 1, j, k, n) * 2.0 + fdata(i - 1, j, k + 1, n)
                            + fdata(i, j, k - 1, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i, j, k + 1, n) * 2.0
                            + fdata(i + 1, j, k - 1, n) + fdata(i + 1, j, k, n) * 2.0 + fdata(i + 1, j, k + 1, n)) / 16.0;
                    else if (K == lo.z || K == hi.z) // Z face
                        cdata(I, J, K, n) =
                        (fdata(i - 1, j - 1, k, n) + fdata(i, j - 1, k, n) * 2.0 + fdata(i + 1, j - 1, k, n)
                            + fdata(i - 1, j, k, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i + 1, j, k, n) * 2.0
                            + fdata(i - 1, j + 1, k, n) + fdata(i, j + 1, k, n) * 2.0 + fdata(i + 1, j + 1, k, n)) / 16.0;
                    else // Interior
                        cdata(I, J, K, n) =
                        (fdata(i - 1, j - 1, k - 1, n) + fdata(i - 1, j - 1, k + 1, n) + fdata(i - 1, j + 1, k - 1, n) + fdata(i - 1, j + 1, k + 1, n) +
                            fdata(i + 1, j - 1, k - 1, n) + fdata(i + 1, j - 1, k + 1, n) + fdata(i + 1, j + 1, k - 1, n) + fdata(i + 1, j + 1, k + 1, n)) / 64.0
                        +
                        (fdata(i, j - 1, k - 1, n) + fdata(i, j - 1, k + 1, n) + fdata(i, j + 1, k - 1, n) + fdata(i, j + 1, k + 1, n) +
                            fdata(i - 1, j, k - 1, n) + fdata(i + 1, j, k - 1, n) + fdata(i - 1, j, k + 1, n) + fdata(i + 1, j, k + 1, n) +
                            fdata(i - 1, j - 1, k, n) + fdata(i - 1, j + 1, k, n) + fdata(i + 1, j - 1, k, n) + fdata(i + 1, j + 1, k, n)) / 32.0
                        +
                        (fdata(i - 1, j, k, n) + fdata(i, j - 1, k, n) + fdata(i, j, k - 1, n) +
                            fdata(i + 1, j, k, n) + fdata(i, j + 1, k, n) + fdata(i, j, k + 1, n)) / 16.0
                        +
                        fdata(i, j, k, n) / 8.0;

#ifdef AMREX_DEBUG
#ifndef ALAMO_GPU
                    if (cdata(I, J, K, n).contains_nan()) Util::Abort(INFO, "restricted model is nan at (", i, ",", j, ",", k, "), fine_amrlev=", fine_amrlev);
#endif
#endif
                }

            });
        }
    }

    crse_ddw.ParallelCopy(fine_ddw_for_coarse, 0, 0, ncomp, 0, 0, cgeom.periodicity());

    FillBoundaryCoeff(crse_ddw,Geom(fine_amrlev,0).periodicity());

}



template<int SYM>
void
Elastic<SYM>::averageDownCoeffsSameAmrLevel(int amrlev)
{
    BL_PROFILE("Elastic::averageDownCoeffsSameAmrLevel()");

    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        amrex::Box cdomain(m_geom[amrlev][mglev].growPeriodicDomain(2));
        cdomain.convert(amrex::IntVect::TheNodeVector());
        amrex::Box fdomain(m_geom[amrlev][mglev - 1].Domain());
        fdomain.convert(amrex::IntVect::TheNodeVector());

        MultiTab& crse = *m_ddw_mf[amrlev][mglev];
        MultiTab& fine = *m_ddw_mf[amrlev][mglev - 1];
        const int ncomp = crse.nComp();

        amrex::BoxArray crseba = crse.boxArray();
        amrex::BoxArray fineba = fine.boxArray();

        BoxArray newba = crseba;
        newba.refine(2);
        MultiTab fine_on_crseba;
        fine_on_crseba.define(newba, crse.DistributionMap(), ncomp, 2);
        fine_on_crseba.ParallelCopy(fine, 0, 0, ncomp, 2, 2,
            m_geom[amrlev][mglev-1].periodicity());
        /* fine_on_crseba.FillBoundaryAndSync(m_geom[amrlev][mglev-1].periodicity()); */

        for (MFIter mfi(crse, false); mfi.isValid(); ++mfi)
        {

            // Restrict valid coarse nodes only; support ghosts extended below after restriction.
            Box bx = mfi.nodaltilebox() & cdomain;

            amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM, SYM>> const& fdata = fine_on_crseba.array(mfi);
            amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& cdata = crse.array(mfi);

            const Dim3 lo = amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

            // (I,J,K) = coarse, (i,j,k) = fine
            amrex::ParallelFor(bx, ncomp, [=] AMREX_GPU_DEVICE(int I, int J, int K, int n) {
                int i = 2 * I, j = 2 * J, k = 2 * K;

                if (n > 0)
                {
                    const int face = n - 1;
                    cdata(I, J, K, n) = 0.5 * (
                        fdata(i, j, k, n)
                        + fdata(i + (face == 0), j + (face == 1),
                            k + (face == 2), n));
                    return;
                }

                if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // Corner
                    cdata(I, J, K, n) = fdata(i, j, k, n);
                else if ((J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // X edge
                    cdata(I, J, K, n) = fdata(i - 1, j, k, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i + 1, j, k, n) * 0.25;
                else if ((K == lo.z || K == hi.z) &&
                    (I == lo.x || I == hi.x)) // Y edge
                    cdata(I, J, K, n) = fdata(i, j - 1, k, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i, j + 1, k, n) * 0.25;
                else if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y)) // Z edge
                    cdata(I, J, K, n) = fdata(i, j, k - 1, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i, j, k + 1, n) * 0.25;
                else if (I == lo.x || I == hi.x) // X face
                    cdata(I, J, K, n) =
                    (fdata(i, j - 1, k - 1, n) + fdata(i, j, k - 1, n) * 2.0 + fdata(i, j + 1, k - 1, n)
                        + fdata(i, j - 1, k, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i, j + 1, k, n) * 2.0
                        + fdata(i, j - 1, k + 1, n) + fdata(i, j, k + 1, n) * 2.0 + fdata(i, j + 1, k + 1, n)) / 16.0;
                else if (J == lo.y || J == hi.y) // Y face
                    cdata(I, J, K, n) =
                    (fdata(i - 1, j, k - 1, n) + fdata(i - 1, j, k, n) * 2.0 + fdata(i - 1, j, k + 1, n)
                        + fdata(i, j, k - 1, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i, j, k + 1, n) * 2.0
                        + fdata(i + 1, j, k - 1, n) + fdata(i + 1, j, k, n) * 2.0 + fdata(i + 1, j, k + 1, n)) / 16.0;
                else if (K == lo.z || K == hi.z) // Z face
                    cdata(I, J, K, n) =
                    (fdata(i - 1, j - 1, k, n) + fdata(i, j - 1, k, n) * 2.0 + fdata(i + 1, j - 1, k, n)
                        + fdata(i - 1, j, k, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i + 1, j, k, n) * 2.0
                        + fdata(i - 1, j + 1, k, n) + fdata(i, j + 1, k, n) * 2.0 + fdata(i + 1, j + 1, k, n)) / 16.0;
                else // Interior
                    cdata(I, J, K, n) =
                    (fdata(i - 1, j - 1, k - 1, n) + fdata(i - 1, j - 1, k + 1, n) + fdata(i - 1, j + 1, k - 1, n) + fdata(i - 1, j + 1, k + 1, n) +
                        fdata(i + 1, j - 1, k - 1, n) + fdata(i + 1, j - 1, k + 1, n) + fdata(i + 1, j + 1, k - 1, n) + fdata(i + 1, j + 1, k + 1, n)) / 64.0
                    +
                    (fdata(i, j - 1, k - 1, n) + fdata(i, j - 1, k + 1, n) + fdata(i, j + 1, k - 1, n) + fdata(i, j + 1, k + 1, n) +
                        fdata(i - 1, j, k - 1, n) + fdata(i + 1, j, k - 1, n) + fdata(i - 1, j, k + 1, n) + fdata(i + 1, j, k + 1, n) +
                        fdata(i - 1, j - 1, k, n) + fdata(i - 1, j + 1, k, n) + fdata(i + 1, j - 1, k, n) + fdata(i + 1, j + 1, k, n)) / 32.0
                    +
                    (fdata(i - 1, j, k, n) + fdata(i, j - 1, k, n) + fdata(i, j, k - 1, n) +
                        fdata(i + 1, j, k, n) + fdata(i, j + 1, k, n) + fdata(i, j, k + 1, n)) / 16.0
                    +
                    fdata(i, j, k, n) / 8.0;

#ifdef AMREX_DEBUG
#ifndef ALAMO_GPU
                if (cdata(I, J, K, n).contains_nan()) Util::Abort(INFO, "restricted model is nan at crse coordinates (I=", I, ",J=", J, ",K=", k, "), amrlev=", amrlev, " interpolating from mglev", mglev - 1, " to ", mglev);
#endif
#endif
            });
        }
        // Extend valid coarse coefficients to smoother support ghosts; FillBoundaryCoeff overwrites with real data.
        for (MFIter mfi(crse, false); mfi.isValid(); ++mfi)
        {
            const Box valid = mfi.validbox();
            const Box grown = mfi.grownnodaltilebox() & cdomain;
            const Dim3 vlo = amrex::lbound(valid), vhi = amrex::ubound(valid);
            amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& cdata = crse.array(mfi);
            amrex::ParallelFor(grown, ncomp,
                [=] AMREX_GPU_DEVICE(int I, int J, int K, int n)
                {
                    if (valid.contains(I, J, K)) return;
                    const int Ic = amrex::max(vlo.x, amrex::min(I, vhi.x));
                    const int Jc = amrex::max(vlo.y, amrex::min(J, vhi.y));
                    const int Kc = amrex::max(vlo.z, amrex::min(K, vhi.z));
                    cdata(I, J, K, n) = cdata(Ic, Jc, Kc, n);
                });
        }
        FillBoundaryCoeff(crse, Geom(amrlev,mglev).periodicity());


        if (!m_psi_set) continue;

        amrex::Box cdomain_cell(m_geom[amrlev][mglev].Domain());
        amrex::Box fdomain_cell(m_geom[amrlev][mglev - 1].Domain());
        MultiFab& crse_psi = *m_psi_mf[amrlev][mglev];
        MultiFab& fine_psi = *m_psi_mf[amrlev][mglev - 1];
        MultiFab fine_psi_on_crseba;
        fine_psi_on_crseba.define(newba.convert(amrex::IntVect::TheCellVector()), crse_psi.DistributionMap(), 1, 1);
        fine_psi_on_crseba.ParallelCopy(fine_psi, 0, 0, 1, 1, 1, m_geom[amrlev][mglev].periodicity());

        for (MFIter mfi(crse_psi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            bx = bx & cdomain_cell;

            amrex::Array4<const Set::Scalar> const& fdata = fine_psi_on_crseba.array(mfi);
            amrex::Array4<Set::Scalar> const& cdata = crse_psi.array(mfi);

            const Dim3 lo = amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

            // (I,J,K) = coarse, (i,j,k) = fine
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int I, int J, int K) {
                int i = 2 * I, j = 2 * J, k = 2 * K;

                if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // Corner
                    cdata(I, J, K) = fdata(i, j, k);
                else if ((J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // X edge
                    cdata(I, J, K) = fdata(i - 1, j, k) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i + 1, j, k) * 0.25;
                else if ((K == lo.z || K == hi.z) &&
                    (I == lo.x || I == hi.x)) // Y edge
                    cdata(I, J, K) = fdata(i, j - 1, k) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i, j + 1, k) * 0.25;
                else if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y)) // Z edge
                    cdata(I, J, K) = fdata(i, j, k - 1) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i, j, k + 1) * 0.25;
                else if (I == lo.x || I == hi.x) // X face
                    cdata(I, J, K) =
                    (fdata(i, j - 1, k - 1) + fdata(i, j, k - 1) * 2.0 + fdata(i, j + 1, k - 1)
                        + fdata(i, j - 1, k) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i, j + 1, k) * 2.0
                        + fdata(i, j - 1, k + 1) + fdata(i, j, k + 1) * 2.0 + fdata(i, j + 1, k + 1)) / 16.0;
                else if (J == lo.y || J == hi.y) // Y face
                    cdata(I, J, K) =
                    (fdata(i - 1, j, k - 1) + fdata(i - 1, j, k) * 2.0 + fdata(i - 1, j, k + 1)
                        + fdata(i, j, k - 1) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i, j, k + 1) * 2.0
                        + fdata(i + 1, j, k - 1) + fdata(i + 1, j, k) * 2.0 + fdata(i + 1, j, k + 1)) / 16.0;
                else if (K == lo.z || K == hi.z) // Z face
                    cdata(I, J, K) =
                    (fdata(i - 1, j - 1, k) + fdata(i, j - 1, k) * 2.0 + fdata(i + 1, j - 1, k)
                        + fdata(i - 1, j, k) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i + 1, j, k) * 2.0
                        + fdata(i - 1, j + 1, k) + fdata(i, j + 1, k) * 2.0 + fdata(i + 1, j + 1, k)) / 16.0;
                else // Interior
                    cdata(I, J, K) =
                    (fdata(i - 1, j - 1, k - 1) + fdata(i - 1, j - 1, k + 1) + fdata(i - 1, j + 1, k - 1) + fdata(i - 1, j + 1, k + 1) +
                        fdata(i + 1, j - 1, k - 1) + fdata(i + 1, j - 1, k + 1) + fdata(i + 1, j + 1, k - 1) + fdata(i + 1, j + 1, k + 1)) / 64.0
                    +
                    (fdata(i, j - 1, k - 1) + fdata(i, j - 1, k + 1) + fdata(i, j + 1, k - 1) + fdata(i, j + 1, k + 1) +
                        fdata(i - 1, j, k - 1) + fdata(i + 1, j, k - 1) + fdata(i - 1, j, k + 1) + fdata(i + 1, j, k + 1) +
                        fdata(i - 1, j - 1, k) + fdata(i - 1, j + 1, k) + fdata(i + 1, j - 1, k) + fdata(i + 1, j + 1, k)) / 32.0
                    +
                    (fdata(i - 1, j, k) + fdata(i, j - 1, k) + fdata(i, j, k - 1) +
                        fdata(i + 1, j, k) + fdata(i, j + 1, k) + fdata(i, j, k + 1)) / 16.0
                    +
                    fdata(i, j, k) / 8.0;
            });
        }
        FillBoundaryCoeff(crse_psi, Geom(amrlev,mglev).periodicity());

    }
}

template<int SYM>
void
Elastic<SYM>::FillBoundaryCoeff(MultiTab& sigma, const amrex::Periodicity& p)
{
    sigma.setMultiGhost(true);
    sigma.FillBoundaryAndSync(p);
}

template<int SYM>
void
Elastic<SYM>::FillBoundaryCoeff(MultiFab& psi, const amrex::Periodicity& p)
{
    psi.setMultiGhost(true);
    psi.FillBoundaryAndSync(p);
}

template class Elastic<Set::Sym::Major>;
template class Elastic<Set::Sym::Isotropic>;
template class Elastic<Set::Sym::MajorMinor>;
template class Elastic<Set::Sym::Diagonal>;

}
