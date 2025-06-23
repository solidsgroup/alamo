#if AMREX_SPACEDIM==2

#ifndef INTEGRATOR_SFI_H
#define INTEGRATOR_SFI_H

#include "Hydro.H"

namespace Integrator
{
template<class PF>
class SFI : 
    virtual public PF,
    virtual public Hydro
{
public:
    SFI():
        PF(),
        Hydro()
    {
    }
    SFI(IO::ParmParse &pp) : SFI()
    {Parse(*this,pp);}
    static void Parse(SFI &value, IO::ParmParse &pp)
    {
        pp.queryclass<PF>(PF::name,&value);
        pp.queryclass<Hydro>("hydro",&value);

        // time to activate hydro integrator
        pp.query_default("tstart",value.tstart,0.0);
        // If true, set hydro_eta to 1-pf_eta
        pp.query_default("invert",value.invert,false);

    }


protected:
    void Initialize(int lev) override
    {
        PF::Initialize(lev);
        Hydro::Initialize(lev);
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        if (a_time > tstart && !hydro_initialized)
        {
            for (int lev = 0; lev <= finest_level; lev++)
            {
                Hydro::Initialize(lev);
                Hydro::Mix(lev);
            }
            hydro_initialized=true;
        }
        PF::TimeStepBegin(a_time, a_step);
        if (hydro_initialized) Hydro::TimeStepBegin(a_time, a_step);
    }

    void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt) override
    {
        PF::Advance(a_lev, a_time, a_dt);
        if (hydro_initialized) Hydro::Advance(a_lev, a_time, a_dt);
    }

    void TagCellsForRefinement(int a_lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        PF::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
        if (hydro_initialized) Hydro::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
    }

    void UpdateEta(int lev, Set::Scalar /*time*/) override
    {
        for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.growntilebox();
            Set::Patch<Set::Scalar> eta_old          = Hydro::eta_old_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> eta_new          = Hydro::eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> pf_eta_old = PF::eta_old_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> pf_eta_new = PF::eta_mf.Patch(lev,mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (invert)
                {
                    eta_new(i,j,k) = 1.0 - pf_eta_new(i,j,k);
                    eta_old(i,j,k) = 1.0 - pf_eta_old(i,j,k);
                }
                else
                {
                    eta_new(i,j,k) = pf_eta_new(i,j,k);
                    eta_old(i,j,k) = pf_eta_old(i,j,k);
                }
            });
        }
        Util::RealFillBoundary(*eta_mf[lev], geom[lev]);
        Util::RealFillBoundary(*eta_old_mf[lev], geom[lev]);
    }

    int invert = false;
    Set::Scalar tstart = NAN;
    bool hydro_initialized = false;
};
} // namespace Integrator
#endif

#endif
