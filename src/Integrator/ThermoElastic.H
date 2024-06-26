#ifndef INTEGRATOR_THERMOELASTIC_H
#define INTEGRATOR_THERMOELASTIC_H

#include "Model/Solid/Affine/Isotropic.H"
#include "Integrator/Mechanics.H"
#include "Integrator/HeatConduction.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
class ThermoElastic : 
    virtual public HeatConduction,
    virtual public Mechanics<Model::Solid::Affine::Isotropic>
{
public:
    ThermoElastic():
        HeatConduction(3),
        Mechanics<Model::Solid::Affine::Isotropic>()
    { }
    ThermoElastic(IO::ParmParse &pp) : ThermoElastic()
    {Parse(*this,pp);}
    static void Parse(ThermoElastic &value, IO::ParmParse &pp)
    {
        pp_queryclass("hc",static_cast<HeatConduction*>(&value));
        pp_queryclass("el",static_cast<Mechanics<Model::Solid::Affine::Isotropic>*>(&value));

        pp_queryarr("alpha",value.alpha); // Diffusion coefficient
    }

protected:
    void Initialize(int lev) override
    {
        HeatConduction::Initialize(lev);
        Mechanics<Model::Solid::Affine::Isotropic>::Initialize(lev);
    }

    void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        Mechanics<Model::Solid::Affine::Isotropic>::UpdateModel(a_step, a_time);

        //Set::Scalar alpha[2];
        //alpha[0] = 0.001; alpha[1] = 0.002;
        
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<Model::Solid::Affine::Isotropic> const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &temp = temp_old_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix F0 = Set::Matrix::Zero();
                    Set::Scalar tempavg = Numeric::Interpolate::CellToNodeAverage(temp,i,j,k,0);
                    for (int n = 0; n < eta.nComp(); n++)
                    {
                        F0 += (eta(i,j,k,n) * alpha[n]) * tempavg * Set::Matrix::Identity();
                    }
                    model(i, j, k).F0 = F0;
                });
            }

            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }

    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        HeatConduction::TimeStepBegin(a_time, a_step);
        Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time, a_step);
    }

    void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt) override
    {
        HeatConduction::Advance(a_lev, a_time, a_dt);
        Mechanics<Model::Solid::Affine::Isotropic>::Advance(a_lev, a_time, a_dt);
    }

    void TagCellsForRefinement(int a_lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        HeatConduction::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
        Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
    }

    std::vector<Set::Scalar> alpha;
};
} // namespace Integrator
#endif
