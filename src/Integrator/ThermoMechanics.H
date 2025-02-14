#ifndef INTEGRATOR_THERMOMECHANICS_H
#define INTEGRATOR_THERMOMECHANICS_H

#include "Model/Solid/Affine/Cubic.H"
#include "Integrator/Mechanics.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Numeric/Stencil.H"

// This is copied from ThermoElastic, adding in cubic model instead of Isotropic.

namespace Integrator
{
class ThermoMechanics :
    virtual public HeatConduction,
    virtual public PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>
{
public:
	ThermoMechanics():
		HeatConduction(3),
		PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>()
	{ }
	ThermoMechanics(IO::ParmParse &pp) : ThermoMechanics()
  {Parse(*this,pp);}
	static void Parse(ThermoMechanics &value, IO::ParmParse &pp)
	{
		pp.queryclass<HeatConduction>("hc",value);
		pp.queryclass<PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>>("pfm",value);
		
		pp_queryarr("alpha",value.alpha); // Diffusion coefficient (heat conduction?)
	}

protected:
    void Initialize(int lev) override
    {
        HeatConduction::Initialize(lev);
        PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>::Initialize(lev);
    }

    void UpdateModel(int a_step, Set::Scalar a_time) override
    {
        PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>::UpdateModel(a_step, a_time);

        //Set::Scalar alpha[2];
        //alpha[0] = 0.001; alpha[1] = 0.002;
        
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                Set::Patch<Model::Solid::Affine::Cubic> &model = model_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> &eta = PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>::eta_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> &temp = HeatConduction::temp_old_mf.Patch(lev,mfi);
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
        PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>::TimeStepBegin(a_time, a_step);
    }

    void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt) override
    {
        HeatConduction::Advance(a_lev, a_time, a_dt);
        PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>::Advance(a_lev, a_time, a_dt);
    }

    void TagCellsForRefinement(int a_lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        HeatConduction::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
        PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>::TagCellsForRefinement(a_lev, a_tags, a_time, a_ngrow);
    }

    std::vector<Set::Scalar> alpha;
};
} // namespace Integrator
#endif
