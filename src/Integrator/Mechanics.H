#ifndef INTEGRATOR_MECHANICS_H
#define INTEGRATOR_MECHANICS_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <cmath>

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
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "IC/PSRead.H"
#include "IC/Constant.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"

#include "Operator/Operator.H"

namespace Integrator
{

template<class MODEL>
class Mechanics : virtual public Base::Mechanics<MODEL>
{
public:
    Mechanics() : Base::Mechanics<MODEL>() { }
    Mechanics(IO::ParmParse& pp) : Base::Mechanics<MODEL>()
    {
        Parse(*this, pp);
    }
    ~Mechanics()
    {
        // Clean up pointers
        delete ic_eta;
        delete ic_psi;
        delete ic_trac_normal;
        delete bc_psi;
        delete bc_trac_normal;
    }

    // ---------------------------------------------------------------------
    // Parse(...) - read parameters from the input
    // ---------------------------------------------------------------------
    static void Parse(Mechanics& value, IO::ParmParse& pp)
    {
        // Let the base class parse its own parameters.
        pp.queryclass<Base::Mechanics<MODEL>>(value);

        // Read number of models.
        int nmodels;
        pp_query_default("nmodels", nmodels, 1);
        for (int i = 0; i < nmodels; i++)
        {
            bool init = false;
            MODEL tmp_model;
            // Try "model"
            if (pp.contains("model"))
            {
                pp_queryclass("model", tmp_model);
                init = true;
            }
            // Try "model1", "model2", etc.
            std::string name = "model" + std::to_string(i + 1);
            if (pp.contains(name.c_str()))
            {
                pp_queryclass(name, tmp_model);
                init = true;
            }
            if (!init)
                Util::Exception(INFO, "No way to initialize ", pp.full(name));
            value.models.push_back(tmp_model);
        }
        Util::Assert(INFO, TEST(value.models.size() > 0));

        // Register "eta" field.
        value.RegisterNodalFab(value.eta_mf, value.models.size(), 2, "eta", true);

        // Read refinement thresholds.
        pp_query_default("eta_ref_threshold", value.m_eta_ref_threshold, 0.01);
        pp_query_default("ref_threshold", value.m_elastic_ref_threshold, 0.01);
        pp_query_default("model_neumann_boundary", value.model_neumann_boundary, false);

        // If multiple models are used and "ic.type" is provided, select the initial condition for eta.
        if (nmodels > 1 && pp.contains("ic.type"))
        {
            pp.select<IC::Constant, IC::Ellipse, IC::Voronoi, IC::BMP, IC::PNG, IC::Expression, IC::PSRead>("ic", value.ic_eta, value.geom);
            value.eta_reset_on_regrid = true;
            pp_query_default("eta.reset_on_regrid", value.eta_reset_on_regrid, true);
        }

        // For psi field.
        if (pp.contains("psi.ic.type"))
        {
            pp.select<IC::Ellipse, IC::Constant, IC::Expression, IC::PSRead, IC::PNG>("psi.ic", value.ic_psi, value.geom);
            value.bc_psi = new BC::Nothing();
            value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", value.plot_psi);
            value.psi_on = true;
            value.psi_reset_on_regrid = true;
            pp_query_default("psi.reset_on_regrid", value.psi_reset_on_regrid, true);
        }

        // For traction-normal field.
        if (pp.contains("trac_normal.ic.type"))
        {
            pp.select<IC::Ellipse, IC::Constant, IC::Expression, IC::PSRead>("trac_normal.ic", value.ic_trac_normal, value.geom);
            if (value.ic_trac_normal)
            {
                value.bc_trac_normal = new BC::Nothing();
                value.RegisterNewFab(value.trac_normal_mf, value.bc_trac_normal, 1, 2, "trac_normal", true);
            }
        }

        // [DAMAGE ADDED] Parse damage parameters.
        pp_query_default("beta_interface", value.beta_interface, 0.0);
        int formula_choice;
        pp_query_default("damage_formula", formula_choice, 2);
        if (formula_choice != 1 && formula_choice != 2)
            formula_choice = 2;
        value.damage_formula = formula_choice;
        pp_query_default("apply_damage_to_stress", value.apply_damage_to_stress, 0);
        value.RegisterNodalFab(value.damage_mf, 1, 2, "damage", true);

        Util::Message(INFO);
    }

    // ---------------------------------------------------------------------
    // Initialize(...) - sets initial conditions.
    // ---------------------------------------------------------------------
    void Initialize(int lev) override
    {
        Base::Mechanics<MODEL>::Initialize(lev);
        eta_mf[lev]->setVal(0.0);
        if (models.size() > 1 && ic_eta)
            ic_eta->Initialize(lev, eta_mf);
        else
            eta_mf[lev]->setVal(1.0);

        if (psi_on)
            ic_psi->Initialize(lev, psi_mf);
        if (ic_trac_normal)
            ic_trac_normal->Initialize(lev, trac_normal_mf);
        damage_mf[lev]->setVal(0.0);
    }

    // ---------------------------------------------------------------------
    // UpdateModel(...) - called each time step to update the effective model.
    // ---------------------------------------------------------------------
    virtual void UpdateModel(int a_step, Set::Scalar time) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable) return;

        if (ic_trac_normal)
        {
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                ic_trac_normal->Initialize(lev, trac_normal_mf, time);
                psi_mf[lev]->FillBoundary();
                amrex::Box domain = geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());
                Set::Vector DX(geom[lev].CellSize());
                for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.nodaltilebox() & domain;
                    auto const& rhs = rhs_mf[lev]->array(mfi);
                    auto const& psi = psi_mf[lev]->array(mfi);
                    auto const& trac_normal = trac_normal_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Vector grad = Numeric::CellGradientOnNode(psi, i, j, k, 0, DX.data());
                        rhs(i,j,k) = trac_normal(i,j,k) * grad;
                    });
                }
                Util::RealFillBoundary(*rhs_mf[lev], geom[lev]);
            }
        }

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            eta_mf[lev]->FillBoundary();
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            Set::Vector DX(geom[lev].CellSize());
            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                auto const& model = model_mf[lev]->array(mfi);
                auto const& eta = eta_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    model(i,j,k) = MODEL::Zero();
                    for (unsigned int n = 0; n < models.size(); n++)
                        model(i,j,k) += eta(i,j,k,n) * models[n];
                });
            }
            if (model_neumann_boundary)
            {
                Util::AssertException(INFO, TEST(AMREX_SPACEDIM == 2),
                                      "neumann boundary works in 2d only");
                for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.grownnodaltilebox() & domain;
                    auto const& model = model_mf[lev]->array(mfi);
                    const Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        if (i == lo.x && j == lo.y)
                            model(i,j,k) = 0.5 * (model(i+1,j,k) + model(i,j+1,k));
                        else if (i == lo.x && j == hi.y)
                            model(i,j,k) = 0.5 * (model(i+1,j,k) + model(i,j-1,k));
                        else if (i == hi.x && j == lo.y)
                            model(i,j,k) = 0.5 * (model(i-1,j,k) + model(i,j+1,k));
                        else if (i == hi.x && j == hi.y)
                            model(i,j,k) = 0.5 * (model(i-1,j,k) + model(i,j-1,k));
                        else if (i == lo.x)
                            model(i,j,k) = model(i+1,j,k);
                        else if (i == hi.x)
                            model(i,j,k) = model(i-1,j,k);
                        else if (j == lo.y)
                            model(i,j,k) = model(i,j+1,k);
                        else if (j == hi.y)
                            model(i,j,k) = model(i,j-1,k);
                    });
                }
            }
            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
        }
    }

    // ---------------------------------------------------------------------
    // Advance(...) - performs the elasticity solve and damage update.
    // ---------------------------------------------------------------------
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override
    {
        const int maxCouplingIter = 3;
        const double tolDamage = 1e-6;
        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());
        double globalMaxDiff = 0.0;
        for (int iter = 0; iter < maxCouplingIter; ++iter)
        {
            // 1) Solve elasticity with the current damage field.
            Base::Mechanics<MODEL>::Advance(lev, time, dt);
            // 2) Damage update.
            double iterMaxDiff = 0.0;
            for (amrex::MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox() & domain;
                auto damage_arr = damage_mf[lev]->array(mfi);
                auto eta_arr = eta_mf[lev]->array(mfi);
                auto stress_arr = this->stress_mf[lev]->array(mfi);
                for (int k = bx.smallEnd()[2]; k <= bx.bigEnd()[2]; ++k)
                {
                    for (int j = bx.smallEnd()[1]; j <= bx.bigEnd()[1]; ++j)
                    {
                        for (int i = bx.smallEnd()[0]; i <= bx.bigEnd()[0]; ++i)
                        {
                            double oldDamage = damage_arr(i,j,k);
                            double e = eta_arr(i,j,k,0);
                            double tol = 1e-6;
                            double beta_val = 0.0;
                            if (e > tol && e < (1.0 - tol))
                                beta_val = beta_interface;
                            double s_xx = stress_arr(i,j,k)(0,0);
                            double s_yy = stress_arr(i,j,k)(1,1);
                            double s_zz = 0.0;
#if AMREX_SPACEDIM==3
                            s_zz = stress_arr(i,j,k)(2,2);
#endif
                            double P_vM = std::sqrt(s_xx*s_xx + s_yy*s_yy + s_zz*s_zz);
                            double dd = 0.0;
                            if (beta_val > 0.0)
                            {
                                if (damage_formula == 1)
                                    dd = beta_val * oldDamage * P_vM;
                                else
                                    dd = beta_val * P_vM * (1.0 - oldDamage);
                            }
                            double newDamage = oldDamage + dd * dt;
                            damage_arr(i,j,k) = newDamage;
                            double diff = std::abs(newDamage - oldDamage);
                            if (diff > iterMaxDiff)
                                iterMaxDiff = diff;
                        }
                    }
                }
                globalMaxDiff = iterMaxDiff;
            }
            if (globalMaxDiff < tolDamage)
                break;
            // 3) Optionally update stress scaling by (1 - d).
            if (apply_damage_to_stress == 1)
            {
                for (amrex::MFIter mfi(*stress_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.nodaltilebox() & domain;
                    auto stress = stress_mf[lev]->array(mfi);
                    auto damage = damage_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        double d_ = damage(i,j,k);
                        double factor = (1.0 - d_);
                        const double eps = 1e-6;
                        if (factor < eps)
                            factor = eps;
                        for (int p = 0; p < AMREX_SPACEDIM; p++)
                            for (int q = 0; q < AMREX_SPACEDIM; q++)
                                stress(i,j,k)(p,q) *= factor;
                    });
                }
            }
        }
    }

    // ---------------------------------------------------------------------
    // TagCellsForRefinement(...) - mark cells for mesh refinement.
    // ---------------------------------------------------------------------
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        if (m_type == Base::Mechanics<MODEL>::Type::Disable)
            return;
        Base::Mechanics<MODEL>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);
        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            auto const& tags = a_tags.array(mfi);
            auto const& eta = eta_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, bx);
                for (int n = 0; n < eta_mf[lev]->nComp(); n++)
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, n, DX.data(), sten);
                    if (grad.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                        tags(i,j,k) = amrex::TagBox::SET;
                }
            });
            if (psi_on)
            {
                auto const& psi = psi_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    Set::Vector gradpsi = Numeric::Gradient(psi, i, j, k, 0, DX.data(), sten);
                    if (gradpsi.lpNorm<2>() * DXnorm > m_eta_ref_threshold)
                        tags(i,j,k) = amrex::TagBox::SET;
                });
            }
        }
    }

    // ---------------------------------------------------------------------
    // Regrid(...) - reinitialize fields when the mesh is adapted.
    // ---------------------------------------------------------------------
    void Regrid(int lev, Set::Scalar time) override
    {
        if (eta_reset_on_regrid && models.size() > 1 && ic_eta)
            ic_eta->Initialize(lev, eta_mf, time);
        if (psi_reset_on_regrid && ic_psi)
            ic_psi->Initialize(lev, psi_mf, time);
    }

protected:
    // Bring in names from the base class.
    using Base::Mechanics<MODEL>::m_type;
    using Base::Mechanics<MODEL>::finest_level;
    using Base::Mechanics<MODEL>::geom;
    using Base::Mechanics<MODEL>::model_mf;
    using Base::Mechanics<MODEL>::rhs_mf;
    using Base::Mechanics<MODEL>::stress_mf;
    using Base::Mechanics<MODEL>::psi_mf;
    using Base::Mechanics<MODEL>::psi_on;

    // The following members are specific to this integrator.
    std::vector<MODEL> models;
    Set::Field<Set::Scalar> eta_mf;
    IC::IC<Set::Scalar>* ic_eta = nullptr;
    bool eta_reset_on_regrid = false;
    IC::IC<Set::Scalar>* ic_psi = nullptr;
    BC::BC<Set::Scalar>* bc_psi = nullptr;
    bool psi_reset_on_regrid = false;
    IC::IC<Set::Scalar>* ic_trac_normal = nullptr;
    BC::BC<Set::Scalar>* bc_trac_normal = nullptr;
    Set::Field<Set::Scalar> trac_normal_mf;
    Set::Scalar m_eta_ref_threshold = 0.01;
    bool model_neumann_boundary = false;
    Set::Scalar m_elastic_ref_threshold = 0.01;

    // Damage-related fields and parameters.
    Set::Field<Set::Scalar> damage_mf;
    Set::Scalar beta_interface = 0.0;
    int damage_formula = 2; // 2: d_new = d_old + dt*(beta*d_old*P_vM); 1: d_new = d_old + dt*(beta*P_vM*(1-d_old))
    int apply_damage_to_stress = 0;
};

} // namespace Integrator

#endif
