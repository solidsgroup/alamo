#ifndef INTEGRATOR_TOPOP_H
#define INTEGRATOR_TOPOP_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Base/Mechanics.H"


#include "IC/IC.H"
#include "BC/BC.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"
#include "BC/Operator/Elastic/Expression.H"

#include "IC/Ellipse.H"
#include "IC/Voronoi.H"
#include "IC/Constant.H"
#include "IC/BMP.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Operator/Operator.H"


namespace Integrator
{
template<class MODEL>
class TopOp: virtual public Base::Mechanics<MODEL>
{
public:

    TopOp(): Base::Mechanics<MODEL>() {}
    TopOp(IO::ParmParse& pp): Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }

    ~TopOp()
    {
        delete ic_psi;
        delete bc_psi;
    }

    // The mechanics integrator manages the solution of an elastic 
    // solve using the MLMG solver. 
    static void Parse(TopOp& value, IO::ParmParse& pp)
    {
        Base::Mechanics<MODEL>::Parse(value, pp);

        pp_queryclass("model", value.model);
        // Read in IC for psi
        if (pp.contains("psi.ic.type"))
        {
            // Initial condition for psi field
            pp.select<IC::Ellipse,IC::Constant>("psi.ic",value.ic_psi,value.geom);

            value.bc_psi = new BC::Constant(1, pp, "psi.bc");
            value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", true);
            value.RegisterNewFab(value.psi_old_mf, value.bc_psi, 1, 2, "psiold", false);
            value.psi_on = true;
        }
        pp_query_default("eta_ref_threshold", value.m_eta_ref_threshold, 0.01); // Refinement threshold based on eta
        pp_query_default("alpha", value.alpha, 1.0); // :math:`\alpha` parameter
        pp_query_default("beta", value.beta, 1.0); // :math:`\beta` parameter
        pp_query_default("gamma", value.gamma, 1.0); // :math:`\gamma` parameter
        pp_queryclass("L", value.L); // Mobility
        if (pp.contains("volume0frac"))
        {
            if (pp.contains("volume0")) Util::Abort(INFO, "Cannot specify volume0frac and volume0");
            Set::Scalar volumefrac;
            pp_query_default("volume0frac", volumefrac, 0.5); // Prescribed volume fraction
            value.volume0 = volumefrac *
                AMREX_D_TERM((value.geom[0].ProbHi()[0] - value.geom[0].ProbLo()[0]),
                    *(value.geom[0].ProbHi()[1] - value.geom[0].ProbLo()[1]),
                    *(value.geom[0].ProbHi()[2] - value.geom[0].ProbLo()[2]));
        }
        else
            pp_query_default("volume0", value.volume0, 0.5); // Prescribed total vlume

        pp_queryclass("lambda", value.lambda); // Lagrange multiplier (can be interplated)


        value.RegisterIntegratedVariable(&value.volume, "volume");
        value.RegisterIntegratedVariable(&value.w_chem_potential, "chem_potential");
        value.RegisterIntegratedVariable(&value.w_bndry, "bndry");
        value.RegisterIntegratedVariable(&value.w_elastic, "elastic");
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        if (psi_on) ic_psi->Initialize(lev, psi_mf);
        if (psi_on) ic_psi->Initialize(lev, psi_old_mf);
    }

    virtual void UpdateModel(int a_step, Set::Scalar /*a_time*/) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        if (a_step > 0) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            model_mf[lev]->setVal(model);
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            Util::RealFillBoundary(*psi_mf[lev], geom[lev]);
            Util::RealFillBoundary(*psi_old_mf[lev], geom[lev]);
        }

    }


    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        BL_PROFILE("TopOp::Advance");
        Base::Mechanics<MODEL>::Advance(lev, time, dt);
        std::swap(psi_old_mf[lev], psi_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

        Set::Scalar Lnow = L(time);
        Set::Scalar lambdaT = lambda(time);

        for (amrex::MFIter mfi(*psi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(1);
            bx = bx & domain;
            amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& psi = (*psi_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& psinew = (*psi_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar driving_force = 0.0;

                driving_force += alpha * 2.0 * psi(i, j, k) * (2.0 * psi(i, j, k) * psi(i, j, k) - 3.0 * psi(i, j, k) + 1.0);
                driving_force += -beta * Numeric::Laplacian(psi, i, j, k, 0, DX);

                Set::Matrix sig_avg = Numeric::Interpolate::NodeToCellAverage(sig, i, j, k, 0);
                Set::Matrix eps_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k, 0);

                driving_force += -gamma * 0.5 * (sig_avg.transpose() * eps_avg).trace() * psi(i, j, k);

                driving_force += lambdaT * (volume - volume0);

                psinew(i, j, k) = psi(i, j, k) - Lnow * dt * driving_force;
                if (psinew(i, j, k) < 0.0) psinew(i, j, k) = 0.0;
                if (psinew(i, j, k) > 1.0) psinew(i, j, k) = 1.0;
            });
        }
    }

    void Integrate(int amrlev, Set::Scalar time, int step,
        const amrex::MFIter& mfi, const amrex::Box& box) override
    {
        BL_PROFILE("TopOp::Integrate");
        Base::Mechanics<MODEL>::Integrate(amrlev, time, step, mfi, box);

        const amrex::Real* DX = geom[amrlev].CellSize();
        Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

        amrex::Array4<amrex::Real> const& psi = (*psi_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Matrix> const& sig = (*stress_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Matrix> const& eps = (*strain_mf[amrlev]).array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            volume += psi(i, j, k, 0) * dv;
            w_chem_potential += alpha * psi(i, j, k, 0) * psi(i, j, k, 0) * (1. - psi(i, j, k, 0) * psi(i, j, k, 0)) * dv;
            w_bndry += beta * 0.5 * Numeric::Gradient(psi, i, j, k, 0, DX).squaredNorm() * dv;

            Set::Matrix sig_avg = Numeric::Interpolate::NodeToCellAverage(sig, i, j, k, 0);
            Set::Matrix eps_avg = Numeric::Interpolate::NodeToCellAverage(eps, i, j, k, 0);
            w_elastic += gamma * 0.5 * (sig_avg.transpose() * eps_avg).trace() * psi(i, j, k) * dv;
        });
    }


    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            if (psi_on)
            {
                amrex::Array4<Set::Scalar> const& psi = psi_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    {
                        Set::Vector gradpsi = Numeric::Gradient(psi, i, j, k, 0, DX.data(), sten);
                        if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                    }
                });
            }
        }
    }


protected:
    MODEL model;
    IC::IC* ic_psi = nullptr;
    BC::BC<Set::Scalar>* bc_psi = nullptr;
    Set::Scalar m_eta_ref_threshold = NAN;
    Set::Field<Set::Scalar> psi_old_mf;


    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::psi_mf;
    using Base::Mechanics<MODEL>::psi_on;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::strain_mf;


    Set::Scalar alpha = NAN;
    Set::Scalar beta = NAN;
    Set::Scalar gamma = NAN;
    //Set::Scalar L = 1.0;
    Numeric::Interpolator::Linear<Set::Scalar> L;
    Set::Scalar volume0 = NAN;
    //Set::Scalar lambda = 1.0;
    Numeric::Interpolator::Linear<Set::Scalar> lambda;

    Set::Scalar volume = NAN;
    Set::Scalar w_chem_potential = NAN;
    Set::Scalar w_bndry = NAN;
    Set::Scalar w_elastic = NAN;
};










} // namespace Integrator
#endif
