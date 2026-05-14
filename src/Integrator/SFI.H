// #if AMREX_SPACEDIM==2

#ifndef INTEGRATOR_SFI_H
#define INTEGRATOR_SFI_H

#include "Hydro.H"

#include <limits>
#include "Integrator/Integrator.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "IC/IC.H"
#include "Set/Set.H"

#include "Base/Mechanics.H"

namespace Integrator
{
template<class PF>
class SFI : 
    virtual public PF,
    virtual public Hydro,
    virtual public Base::Mechanics<Model::Solid::Finite::NeoHookeanPredeformed>
{
    using model_type = Model::Solid::Finite::NeoHookeanPredeformed;
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

        // Whether to use Neo-hookean Elastic model
        pp.query_default("elastic.on", value.elastic.on, 0); 

        // Body force
        pp.query_default("elastic.traction", value.elastic.traction, 0.0);

        pp.queryclass<Base::Mechanics<model_type>>("elastic",value);
        amrex::Print() << "After Mechanics parse: plot_stress=" << value.plot_stress 
               << " m_type=" << (int)value.m_type << "\n";

        if (value.m_type != Type::Disable)
        {
            // elastic model of Steel
            pp.queryclass<Model::Solid::Finite::NeoHookeanPredeformed>("model_steel", value.elastic.model_steel);
            // Use Hydro's eta field as the psi field for the solver
            value.psi_on = false;
        }

    }


protected:
    void Initialize(int lev) override
    {
        PF::Initialize(lev);
        Hydro::Initialize(lev);
        psi_on = false;
        solver.setPsi(*Hydro::eta_mf);

        Base::Mechanics<model_type>::Initialize(lev);
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {   
        if (a_time > tstart && !hydro_initialized)
        {
            for (int lev = 0; lev <= finest_level; lev++)
            {
                Hydro::Initialize(lev);
                Hydro::Mix(lev);
                Base::Mechanics<model_type>::Initialize(lev);
            }
            hydro_initialized = true;
        }

        PF::TimeStepBegin(a_time, a_step);

        if (hydro_initialized)
        {
            Hydro::TimeStepBegin(a_time, a_step);

            // Refresh psi pointer and fill ghost cells before mechanics solve
            solver.setPsi(*Hydro::eta_mf);
            for (int lev = 0; lev <= finest_level; lev++) {
                // amrex::Print() << "eta nGrow = " << (*Hydro::eta_mf)[lev]->nGrow() << "\n";
                (*Hydro::eta_mf)[lev]->FillBoundary(geom[lev].periodicity());
                // Try this instead - consistent with how model_mf is filled elsewhere in your code
                // Util::RealFillBoundary(*Hydro::eta_mf[lev], geom[lev]);
            }
            Base::Mechanics<model_type>::TimeStepBegin(a_time, a_step);
        }
    }

    void TimeStepComplete(Set::Scalar a_time, int a_iter) override
    {
        PF::TimeStepComplete(a_time, a_iter);
        if (hydro_initialized) Hydro::TimeStepComplete(a_time, a_iter);
    }

    void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt) override
    {
        PF::Advance(a_lev, a_time, a_dt);
        if (hydro_initialized) {
            Hydro::Advance(a_lev, a_time, a_dt);
        }
        Base::Mechanics<model_type>::Advance(a_lev, a_time, a_dt);
    }

    void TagCellsForRefinement(int a_lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        PF::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
        if (hydro_initialized) Hydro::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        PF::Regrid(lev, time);
        if (hydro_initialized) Hydro::Regrid(lev, time);
    }

    void UpdateEta(int lev, Set::Scalar /*time*/) override
    {
        for (amrex::MFIter mfi(*(*eta_mf)[lev], true); mfi.isValid(); ++mfi)
        {
            //const amrex::Box& bx = mfi.growntilebox(); // TODO ask Dr. R why I need this change
            const amrex::Box& bx = mfi.tilebox();
            
            Set::Patch<Set::Scalar> eta_old          = Hydro::eta_old_mf->Patch(lev,mfi);
            Set::Patch<Set::Scalar> eta_new          = Hydro::eta_mf->Patch(lev,mfi);
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
        Util::RealFillBoundary(*(*eta_mf)[lev], geom[lev]);
        Util::RealFillBoundary(*(*eta_old_mf)[lev], geom[lev]);
    }

    void UpdateModel(int /*a_step*/, Set::Scalar /*a_time*/) override
    {   
        // amrex::Print() << "SFI::UpdateModel called, m_type=" << (int)m_type 
        //         << " Type::Static=" << (int)Base::Mechanics<model_type>::Type::Static
        //         << " Type::Disable=" << (int)Base::Mechanics<model_type>::Type::Disable << "\n";
        // amrex::Print() << "model_steel.mu = " << elastic.model_steel.mu << "\n";
        // amrex::Print() << "model_steel.F0 = \n" << elastic.model_steel.F0 << "\n";
        // amrex::Print() << "elastic.traction = " << elastic.traction << "\n";

        if (m_type == Base::Mechanics<model_type>::Type::Disable) 
        {
            amrex::Print() << "RETURNING EARLY\n";
            return;
        }
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            amrex::Box domain = this->geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            const Set::Scalar* DX = geom[lev].CellSize();

            // Match Flame's pattern - FillBoundary on inputs before MFIter
            (*Hydro::eta_mf)[lev]->FillBoundary(geom[lev].periodicity());

            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box smallbox = mfi.nodaltilebox();
                amrex::Box bx = mfi.grownnodaltilebox() & domain; // TODO ask Dr.R about this change

                // Match Flame's Patch() pattern exactly
                Set::Patch<model_type>        model = model_mf.Patch(lev, mfi);
                Set::Patch<Set::Vector>       rhs   = rhs_mf.Patch(lev, mfi);
                Set::Patch<const Set::Scalar> eta   = Hydro::eta_mf->Patch(lev, mfi);  // TODO ask Dr. R about this change
                Set::Patch<const Set::Scalar> pres  = Hydro::pressure_mf.Patch(lev, mfi);

                if (elastic.on)
                {
                    amrex::ParallelFor(smallbox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Vector grad_eta = Numeric::CellGradientOnNode(eta, i, j, k, 0, DX);
                        Set::Scalar pres_avg = Numeric::Interpolate::CellToNodeAverage(pres, i, j, k, 0);
                        rhs(i, j, k) = (elastic.traction) * grad_eta + pres_avg*grad_eta;
                    });
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        model_type model_steel = elastic.model_steel;
                        model_steel.F0 -= Set::Matrix::Identity();
                        model_steel.F0 += Set::Matrix::Identity();
                        model(i, j, k) = model_steel;
                    });
                }
                else
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        model_type model_steel = elastic.model_steel;
                        model_steel.F0 *= Set::Matrix::Zero();
                        model(i, j, k) = model_steel;
                    });
                }
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

    int invert = false;
    Set::Scalar tstart = NAN;
    bool hydro_initialized = false;

    // Thermoelastic variables
    struct {
        int on = 0;
        Set::Scalar Telastic = NAN;
        model_type model_steel;
        Set::Scalar traction = NAN;
    } elastic;
};
} // namespace Integrator
#endif

// #endif
