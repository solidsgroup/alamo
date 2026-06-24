#ifndef OPERATOR_SPECTRAL_FFT_H_
#define OPERATOR_SPECTRAL_FFT_H_

#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <memory>

#include "BC/Nothing.H"
#include "Set/Set.H"
#include "Util/Util.H"

#ifdef ALAMO_FFT
#include <AMReX_FFT.H>
#endif

namespace Operator
{
namespace Spectral
{

#ifdef ALAMO_FFT
class FFT
{
public:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar K (int dir, int mode, const int* length, const Set::Scalar* kscale)
    {
        if (dir == 0) return mode * kscale[0];
        return (mode < length[dir] / 2 ? mode : mode - length[dir]) * kscale[dir];
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar Omega2 (int m, int n, int p, const int* length, const Set::Scalar* kscale)
    {
#if AMREX_SPACEDIM < 3
        amrex::ignore_unused(p);
#endif
        return AMREX_D_TERM(K(0,m,length,kscale) * K(0,m,length,kscale),
            + K(1,n,length,kscale) * K(1,n,length,kscale),
            + K(2,p,length,kscale) * K(2,p,length,kscale));
    }

    FFT (const amrex::Geometry& a_geom)
        : m_geom(a_geom),
          m_fft(a_geom.Domain())
    {
        auto const& layout = m_fft.getSpectralDataLayout();
        m_cba = layout.first;
        m_cdm = layout.second;

        m_domain = a_geom.Domain();
        m_scaling = 1.0 / a_geom.Domain().d_numPts();
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            m_length[dir] = a_geom.Domain().length(dir);
            m_kscale[dir] = 2.0 * Set::Constant::Pi / (m_length[dir] * a_geom.CellSize(dir));
        }
    }

    FFT (
        const amrex::Vector<amrex::Geometry>& a_geom,
        const amrex::Vector<amrex::IntVect>& a_ref_ratio,
        int a_finest_level)
        : FFT(a_geom[a_finest_level])
    {
        Util::Assert(INFO, TEST(a_finest_level >= 0));
        Util::Assert(INFO, TEST(static_cast<int>(a_geom.size()) > a_finest_level));
        Util::Assert(INFO, TEST(a_finest_level == 0 || static_cast<int>(a_ref_ratio.size()) >= a_finest_level));
        m_geom_hierarchy = a_geom;
        m_ref_ratio = a_ref_ratio;
    }

    FFT (const FFT&) = delete;
    FFT (FFT&&) = delete;
    FFT& operator= (const FFT&) = delete;
    FFT& operator= (FFT&&) = delete;

    const amrex::Geometry& Geom () const { return m_geom; }
    const amrex::Box& Domain () const { return m_domain; }
    const amrex::BoxArray& SpectralBoxArray () const { return m_cba; }
    const amrex::DistributionMapping& SpectralDistributionMap () const { return m_cdm; }
    int Length (int dir) const { return m_length[dir]; }
    Set::Scalar KScale (int dir) const { return m_kscale[dir]; }
    Set::Scalar K (int dir, int mode) const
    {
        return K(dir, mode, m_length, m_kscale);
    }

    Set::Scalar Omega2 (int m, int n, int p) const
    {
        return Omega2(m, n, p, m_length, m_kscale);
    }

    amrex::FabArray<amrex::BaseFab<Set::Complex> > MakeSpectralFab (int ncomp = 1, int ngrow = 0) const
    {
        return amrex::FabArray<amrex::BaseFab<Set::Complex> >(m_cba, m_cdm, ncomp, ngrow);
    }

    template <class F>
    void ParallelFor (const amrex::Box& bx, F&& f) const
    {
        int length[AMREX_SPACEDIM];
        Set::Scalar kscale[AMREX_SPACEDIM];
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            length[dir] = m_length[dir];
            kscale[dir] = m_kscale[dir];
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {
            f(m, n, p, Omega2(m, n, p, length, kscale));
        });
    }

    template <class F>
    void ParallelForModes (const amrex::Box& bx, F&& f) const
    {
        int length[AMREX_SPACEDIM];
        Set::Scalar kscale[AMREX_SPACEDIM];
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            length[dir] = m_length[dir];
            kscale[dir] = m_kscale[dir];
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int m, int n, int p) {
            Set::Scalar k1 = K(0, m, length, kscale);
            Set::Scalar k2 = 0.0;
            Set::Scalar k3 = 0.0;
#if AMREX_SPACEDIM >= 2
            k2 = K(1, n, length, kscale);
#endif
#if AMREX_SPACEDIM >= 3
            k3 = K(2, p, length, kscale);
#endif
            Set::Scalar omega2 = k1*k1 + k2*k2 + k3*k3;
            f(m, n, p, k1, k2, k3, omega2);
        });
    }

    void Forward (const amrex::MultiFab& in, amrex::FabArray<amrex::BaseFab<Set::Complex> >& out, int scomp = 0, int dcomp = 0)
    {
        m_fft.forward(in, out, scomp, dcomp);
    }

    void Backward (const amrex::FabArray<amrex::BaseFab<Set::Complex> >& in, amrex::MultiFab& out, int scomp = 0, int dcomp = 0)
    {
        m_fft.backward(in, out, scomp, dcomp);
        out.mult(m_scaling, dcomp, 1, 0);
    }

    bool IsFullDomain (const amrex::MultiFab& mf) const
    {
        return mf.boxArray().complementIn(m_domain).isEmpty();
    }

    void Forward (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& hierarchy,
        int lev,
        amrex::FabArray<amrex::BaseFab<Set::Complex> >& out,
        int scomp = 0,
        int dcomp = 0,
        amrex::Real time = 0.0)
    {
        Util::Assert(INFO, TEST(lev >= 0));
        Util::Assert(INFO, TEST(static_cast<int>(hierarchy.size()) > lev));
        if (IsFullDomain(*hierarchy[lev]))
        {
            Forward(*hierarchy[lev], out, scomp, dcomp);
            return;
        }

        auto full = CompositeToUniform(hierarchy, lev, time, hierarchy[lev]->nComp(), hierarchy[lev]->nGrow());
        Forward(*full, out, scomp, dcomp);
    }

    void Backward (
        const amrex::FabArray<amrex::BaseFab<Set::Complex> >& in,
        amrex::Vector<std::unique_ptr<amrex::MultiFab> >& hierarchy,
        int lev,
        int scomp = 0,
        int dcomp = 0)
    {
        Util::Assert(INFO, TEST(lev >= 0));
        Util::Assert(INFO, TEST(static_cast<int>(hierarchy.size()) > lev));
        if (IsFullDomain(*hierarchy[lev]))
        {
            Backward(in, *hierarchy[lev], scomp, dcomp);
            return;
        }

        EnsureRealTemporary(hierarchy[lev]->nComp(), hierarchy[lev]->nGrow());
        Backward(in, *m_real_tmp, scomp, dcomp);
        CopyFromFull(*m_real_tmp, *hierarchy[lev], m_geom, dcomp, dcomp, 1);
    }

    static amrex::BoxArray FullBoxArray (const amrex::Geometry& geom, int max_grid_size = 0)
    {
        amrex::BoxArray ba(geom.Domain());
        if (max_grid_size > 0) ba.maxSize(max_grid_size);
        return ba;
    }

    static std::unique_ptr<amrex::MultiFab> MakeFullMultiFab (
        const amrex::Geometry& geom,
        int ncomp,
        int ngrow = 0,
        int max_grid_size = 0)
    {
        amrex::BoxArray ba = FullBoxArray(geom, max_grid_size);
        amrex::DistributionMapping dm(ba);
        return std::make_unique<amrex::MultiFab>(ba, dm, ncomp, ngrow);
    }

    static void CopyToFull (
        amrex::MultiFab& full,
        const amrex::MultiFab& src,
        const amrex::Geometry& geom,
        int scomp = 0,
        int dcomp = 0,
        int ncomp = -1)
    {
        if (ncomp < 0) ncomp = std::min(src.nComp() - scomp, full.nComp() - dcomp);
        full.ParallelCopy(src, scomp, dcomp, ncomp, 0, 0, geom.periodicity());
    }

    static void CopyFromFull (
        const amrex::MultiFab& full,
        amrex::MultiFab& dst,
        const amrex::Geometry& geom,
        int scomp = 0,
        int dcomp = 0,
        int ncomp = -1)
    {
        if (ncomp < 0) ncomp = std::min(full.nComp() - scomp, dst.nComp() - dcomp);
        dst.ParallelCopy(full, scomp, dcomp, ncomp, 0, 0, geom.periodicity());
    }

    template <class PhysBC>
    static void ProlongateFromCoarse (
        amrex::MultiFab& fine,
        amrex::Real time,
        const amrex::MultiFab& coarse,
        const amrex::Geometry& crse_geom,
        const amrex::Geometry& fine_geom,
        const amrex::IntVect& ref_ratio,
        PhysBC& physbc,
        int scomp = 0,
        int dcomp = 0,
        int ncomp = -1)
    {
        if (ncomp < 0) ncomp = std::min(coarse.nComp() - scomp, fine.nComp() - dcomp);

        amrex::Interpolater* mapper = &amrex::cell_cons_interp;
        if (fine.boxArray().ixType() == amrex::IndexType::TheNodeType())
            mapper = &amrex::node_bilinear_interp;

        physbc.define(fine_geom);
        amrex::Vector<amrex::BCRec> bcs(ncomp, physbc.GetBCRec());
        amrex::InterpFromCoarseLevel(
            fine, time, coarse, scomp, dcomp, ncomp,
            crse_geom, fine_geom,
            physbc, 0,
            physbc, 0,
            ref_ratio,
            mapper, bcs, 0);
    }

    template <class PhysBC>
    static std::unique_ptr<amrex::MultiFab> CompositeToUniform (
        const amrex::Vector<amrex::MultiFab*>& hierarchy,
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::IntVect>& ref_ratio,
        int finest_level,
        amrex::Real time,
        PhysBC& physbc,
        int ncomp,
        int ngrow = 0,
        int max_grid_size = 0)
    {
        Util::Assert(INFO, TEST(finest_level >= 0));
        Util::Assert(INFO, TEST(static_cast<int>(hierarchy.size()) > finest_level));
        Util::Assert(INFO, TEST(static_cast<int>(geom.size()) > finest_level));

        amrex::Vector<std::unique_ptr<amrex::MultiFab> > full(finest_level + 1);
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            full[lev] = MakeFullMultiFab(geom[lev], ncomp, ngrow, max_grid_size);
            full[lev]->setVal(0.0);

            if (lev > 0)
            {
                Util::Assert(INFO, TEST(static_cast<int>(ref_ratio.size()) >= lev));
                ProlongateFromCoarse(*full[lev], time, *full[lev - 1], geom[lev - 1], geom[lev],
                    ref_ratio[lev - 1], physbc, 0, 0, ncomp);
            }

            CopyToFull(*full[lev], *hierarchy[lev], geom[lev], 0, 0, ncomp);
        }

        return std::move(full[finest_level]);
    }

    static void UniformToComposite (
        const amrex::MultiFab& full_fine,
        amrex::Vector<amrex::MultiFab*>& hierarchy,
        const amrex::Vector<amrex::Geometry>& geom,
        const amrex::Vector<amrex::IntVect>& ref_ratio,
        int finest_level,
        int ncomp,
        int max_grid_size = 0)
    {
        Util::Assert(INFO, TEST(finest_level >= 0));
        Util::Assert(INFO, TEST(static_cast<int>(hierarchy.size()) > finest_level));
        Util::Assert(INFO, TEST(static_cast<int>(geom.size()) > finest_level));

        amrex::Vector<std::unique_ptr<amrex::MultiFab> > full(finest_level + 1);
        full[finest_level] = MakeFullMultiFab(geom[finest_level], ncomp, 0, max_grid_size);
        CopyToFull(*full[finest_level], full_fine, geom[finest_level], 0, 0, ncomp);

        for (int lev = finest_level - 1; lev >= 0; --lev)
        {
            Util::Assert(INFO, TEST(static_cast<int>(ref_ratio.size()) > lev));
            full[lev] = MakeFullMultiFab(geom[lev], ncomp, 0, max_grid_size);
            amrex::average_down(*full[lev + 1], *full[lev], 0, ncomp, ref_ratio[lev]);
        }

        for (int lev = 0; lev <= finest_level; ++lev)
            CopyFromFull(*full[lev], *hierarchy[lev], geom[lev], 0, 0, ncomp);
    }

private:
    void EnsureRealTemporary (int ncomp, int ngrow)
    {
        if (!m_real_tmp || m_real_tmp->nComp() != ncomp || m_real_tmp->nGrow() != ngrow)
            m_real_tmp = MakeFullMultiFab(m_geom, ncomp, ngrow);
    }

    std::unique_ptr<amrex::MultiFab> CompositeToUniform (
        const amrex::Vector<std::unique_ptr<amrex::MultiFab> >& hierarchy,
        int finest_lev,
        amrex::Real time,
        int ncomp,
        int ngrow)
    {
        Util::Assert(INFO, TEST(finest_lev >= 0));
        Util::Assert(INFO, TEST(static_cast<int>(hierarchy.size()) > finest_lev));
        Util::Assert(INFO, TEST(static_cast<int>(m_geom_hierarchy.size()) > finest_lev),
            "Patchy spectral FFT requires hierarchy geometry context");
        Util::Assert(INFO, TEST(finest_lev == 0 || static_cast<int>(m_ref_ratio.size()) >= finest_lev));

        amrex::Vector<std::unique_ptr<amrex::MultiFab> > full(finest_lev + 1);
        for (int ilev = 0; ilev <= finest_lev; ++ilev)
        {
            const amrex::Geometry& level_geom = m_geom_hierarchy[ilev];
            full[ilev] = MakeFullMultiFab(level_geom, ncomp, ngrow);
            full[ilev]->setVal(0.0);

            if (ilev > 0)
            {
                BC::Nothing physbc;
                ProlongateFromCoarse(*full[ilev], time, *full[ilev - 1],
                    m_geom_hierarchy[ilev - 1], level_geom, m_ref_ratio[ilev - 1], physbc, 0, 0, ncomp);
            }

            CopyToFull(*full[ilev], *hierarchy[ilev], level_geom, 0, 0, ncomp);
        }

        return std::move(full[finest_lev]);
    }

    amrex::Geometry m_geom;
    amrex::FFT::R2C<Set::Scalar> m_fft;
    amrex::BoxArray m_cba;
    amrex::DistributionMapping m_cdm;
    std::unique_ptr<amrex::MultiFab> m_real_tmp;
    amrex::Vector<amrex::Geometry> m_geom_hierarchy;
    amrex::Vector<amrex::IntVect> m_ref_ratio;
    amrex::Box m_domain;
    int m_length[AMREX_SPACEDIM];
    Set::Scalar m_kscale[AMREX_SPACEDIM];
    Set::Scalar m_scaling;
};
#else
class FFT
{
public:
    FFT (const amrex::Geometry&) { Util::Abort(INFO,"Alamo must be compiled with fft"); }
};
#endif

}
}

#endif
