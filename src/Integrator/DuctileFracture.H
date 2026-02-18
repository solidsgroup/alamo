//
// This class implements second and fourth order phase field brittle fracture with near singular solver.
// This class inherits from the base class :code:`Base/Mechanics.H`.
// See `this link <https://doi.org/10.1007/s00466-023-02325-8>`_
// for more information on this implementation.
//
// The energy functional for a second order model is given by
// .. math::
//
//      \mathcal{L}= \int_\Omega \left[\left(g(c) + \eta \right)W_0^+ + W_0^-\right] dV + \int_\Omega G_c \left[ \frac{w(c)}{4\xi} + \xi |\nabla c|^2\right] dV
//
// where :math:`g(c)` and :math:`w(c)` are interpolation functions specified by the user.
// The fracture energy :math:`G_c` and crack length scale :math:`xi` are also read from input file.
// The tension compression asymmetry is accounted by splitting the strain energy :math:`W_0` as
// .. math::
//
// W_0^\pm = \frac{1}{2} \lambda (\operatorname{tr}\bm{\varepsilon}_\pm)^2 + \mu \operatorname{tr}\left(\bm{\varepsilon}_\pm^2\right), \quad \bm{\varepsilon_\pm} = \sum_{i=1}^d \left(\varepsilon_i \right)_\pm \hat{\bm{v}}_i\otimes\hat{\bm{v}}_i
//
// where :math:`\lambda` and :math:`\mu` are material models for linear elastic isotropic material, :math:`\bm{\varepsilon}_\pm` are the positive and negative components of the strain tensor :math:`\bm{\varepsilon}` computed through eigenvalue decomposition.
// The fourth order model adds a laplacian term to the free energy functional.
// The code performs a staggered solve where it solves the elastic problem implicitly and crack problem explicitly.
//
// Class methods:
//
// #. :code:`Fracture()`:
//    Basic constructor. Does nothing, and leaves all values initiated as NAN.
// #. :code:`Fracture(IO::ParmParse &pp)`:
//    Calls the parser.
// #. :code:`static void Parse(Fracture &value, IO::ParmParse &pp)`
//    Parses input file, ICs, BCs, and sets up multifabs appropriately
// #. :code:`void Initialize(int lev) override`
//    Calls IC and sets up the initial crack geometry.
// #. :code:`virtual void UpdateModel(int a_step) override`
//    Performs degradation by updating the :math:`\psi` field using :math:`g(c)`.
// #. :code:`void TimeStepBegin(Set::Scalar a_time, int a_step) override`
//    Solves the elastic problem, performs eigen value decomposition of strain, and computes the positive part of the strain energy.
// #. :code:`void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt)`
//    Advances the crack field by computing the variational derivative of energy functional with :math:`c`.
// #. :code:`void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar a_time, int a_ngrow) override`
//    Refines the grid based on the crack field.
// #. :code:`void Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &a_box) override`
//    Performs spatial integration of the driving force to check for convergence of crack problem
// #. :code:`void TimeStepComplete(Set::Scalar /*time*/, int /* iter*/)`
//    Checks whether the solver should work on crack problem or elastic problem.
//
#ifndef INTEGRATOR_DUCTILEFRACTURE_H
#define INTEGRATOR_DUCTILEFRACTURE_H

// AMReX Includes
#include "AMReX_Array4.H"
#include "AMReX_GpuComplex.H"
#include "AMReX_MFIter.H"
#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TimeIntegrator.H>

#include "Integrator/Base/Mechanics.H"

#include "Model/Interface/Crack/PFCZM.H"
#include "Model/Solid/Affine/J2Plastic.H"

#include "BC/Constant.H"
#include "IC/BMP.H"
#include "IC/Ellipse.H"
// #include "IC/Ellipsoid.H"
#include "IC/Expression.H"
#include "IC/IC.H"
#include "IC/Laminate.H"
#include "IC/Notch.H"
#include "IC/PNG.H"
#include "IC/PerturbedInterface.H"

#include "Numeric/Stencil.H"

#include <cmath>
#include <eigen3/Eigen/Dense>

namespace Integrator
{
using ductile_model = Model::Solid::Affine::J2Plastic;
using crack_type = Model::Interface::Crack::PFCZM;

class DuctileFracture : virtual public Base::Mechanics<ductile_model>
{
public:
    static constexpr const char *name = "fracture";
    DuctileFracture() : Base::Mechanics<ductile_model>() {}

    DuctileFracture(IO::ParmParse &pp) : Base::Mechanics<ductile_model>()
    {
        Parse(*this, pp);
    }

    static void
    Parse(DuctileFracture &value, IO::ParmParse &pp)
    {
        Base::Mechanics<ductile_model>::Parse(value, pp);

        // Material field related parsing
        pp.query("material.refinement_threshold", value.material.m_eta_ref_threshold);

        // Driving force threshold
        pp.query_default("driving_force_refinement_threshold", value.crack.driving_force_refinement_threshold, 1E100);

        // Let's figure out if we are doing multi-material simulations.
        std::string mat_ic_type;
        if (pp.contains("material.ic.type"))
        {
            // IC determining material distribution
            pp.select<IC::Laminate, IC::Ellipse, IC::PerturbedInterface, IC::BMP, IC::PNG, IC::Expression>("material.ic", value.material.ic, value.geom);
            value.material.is_ic = true;
        }
        else
        {
            value.material.is_ic = false;
        }

        // Let's query material properties for different materials in the system
        pp.queryclass_enumerate<ductile_model>("material.model", value.material.models);
        value.material.num_mat = value.material.models.size();
        Util::Assert(INFO, TEST(value.material.models.size() > 0));

        // This is the fab that stores the material field.
        // For a simple, single material simulation, this will be a uniform field with value 1 everywhere
        // plotting parameters for eta field
        pp.query_default("material.plot_eta", value.material.plot_eta, false);
        value.RegisterNewFab(value.material.eta_mf, &value.bc_nothing, value.material.num_mat, 3, "eta", value.material.plot_eta);

        // Crack related parsing
        pp.query_default("crack.refinement_threshold", value.crack.refinement_threshold, 0.0001); // mesh refinement criteria
        pp.query_default("crack.df.beta", value.crack.beta, 0.0);                                 // constant multiplier for fourth orther model with bilaplacian
        pp.query_default("crack.df.tol_rel", value.crack.tol_rel, 1E-3);                          // relative tolerance for convergence of driving force
        pp.query_default("crack.df.tol_abs", value.crack.tol_abs, 1E-3);                          // absolute tolerance for convergence of driving force
        pp.query_default("crack.df.max_iter", value.crack.max_iter, 500.);                        // Maximum number of solver iterations
        pp.query_default("crack.df.mult_Gc", value.crack.mult_Gc, 1.0);                           // constant multiplier for controlling fracture energy of interface
        pp.query_default("crack.df.mult_lap", value.crack.mult_lap, 1.0);                         // constant multiplier for second order model with laplacian
        pp.query_default("crack.df.el_mult", value.crack.el_mult, 1.0);                           // unit conversion multiplier between elastic problem and crack problem
        pp.query_default("crack.plot_df", value.crack.plot_df, false);
        pp.query_default("crack.plot_energy", value.crack.plot_energy, false);
        pp.query_default("crack.plot_crack_old", value.crack.plot_crack_old, false);

        // Let's figure out if there is an initial crack or void
        pp.select_enumerate<IC::Constant,IC::Notch, IC::Ellipse, IC::Expression>("crack.ic", value.crack.ic, value.geom);

        // value.RegisterNodalFab(value.crack.c_mf, 1, 2, "crack", true);
        // value.RegisterNodalFab(value.crack.c_old_mf, 1, 2, "crack_old", true);
        // value.RegisterNodalFab(value.crack.driving_force_mf, 5, 2, "driving_force", true);
        value.RegisterNodalFab(value.crack.energy_pristine_mf, 3, 2, "energy_pristine", value.crack.plot_energy);
        value.RegisterNodalFab(value.crack.history_var_mf, 3, 2, "history_variable", value.crack.plot_energy);
        value.RegisterIntegratedVariable(&(value.crack.driving_force_norm), "driving_force_norm");
        value.RegisterIntegratedVariable(&(value.crack.crack_l2_err), "crack_l2_err");
        value.RegisterIntegratedVariable(&(value.crack.crack_norm), "crack_norm");

        // By default the psi variable will be on for fracture simulations
        value.psi_on = true;
        value.bc_psi = new BC::Constant(1, pp, "crack.bc"); // boundary conditions for crack
        value.bc_df = new BC::Constant(5, pp, "crack.bc"); // boundary conditions for crack
        
        value.RegisterNewFab(value.crack.c_mf, value.bc_psi, 1, 2, "crack", true);
        value.RegisterNewFab(value.crack.c_old_mf, value.bc_psi, 1, 2, "crack_old", value.crack.plot_crack_old);
        value.RegisterNewFab(value.crack.driving_force_mf,value.bc_df, 5, 2, "driving_force", value.crack.plot_df);
        value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", true);

        // This is needed for the specific crack model we are implementing
        pp.queryclass_enumerate<crack_type>("crack.model",value.crack.cracktype);
        Util::Assert(INFO, TEST(value.material.models.size() > 0));
    }

    void
    Initialize(int lev) override
    {
        Base::Mechanics<ductile_model>::Initialize(lev);

        crack.driving_force_mf[lev]->setVal(0.0);
        crack.c_mf[lev]->setVal(1.0);
        crack.c_old_mf[lev]->setVal(1.0);

        // For now, we are setting psi field to 1. This prevents nans.
        // We will modify the psi field later.
        psi_mf[lev]->setVal(1.0);

        if (crack.is_ic)
        {
            for (unsigned int i = 0; i < crack.ic.size(); i++)
            {
                // This initializes the crack field.
                crack.ic[i]->Add(lev, crack.c_mf);
                crack.ic[i]->Add(lev, crack.c_old_mf);
            }
        }

        // This initializes the material field
        material.eta_mf[lev]->setVal(1.0);
        if (material.is_ic)
            material.ic->Initialize(lev, material.eta_mf);
        else
            material.eta_mf[lev]->setVal(1.0);

        // set history variable to zero
        crack.energy_pristine_mf[lev]->setVal(0.0);
        crack.history_var_mf[lev]->setVal(0.0);
    }

    virtual void
    UpdateModel(int a_step, Set::Scalar /*a_time*/) override
    {
        if (m_type == Base::Mechanics<ductile_model>::Type::Disable)
            return;

        // This sets the model field for the very first time step
        // We should not be resetting model field everytime.
        if (a_step == 0)
        {
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                material.eta_mf[lev]->FillBoundary();

                for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.grownnodaltilebox();
                    amrex::Array4<ductile_model> const &model = model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        model(i, j, k) = ductile_model::Zero();
                        for (int n = 0; n < material.num_mat; n++){
                            // model(i, j, k) += eta(i, j, k, n) * material.models[n];
                            model(i, j, k) += Numeric::Interpolate::CellToNodeAverage(eta, i, j, k, n) * material.models[n];
                        }
                    });
                }
                Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            }
        }

        // This is where we perform "degradation"
        // Essentially we will set the psi field based on the c^2 value.
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            crack.c_mf[lev]->FillBoundary();
            for (MFIter mfi(*psi_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.growntilebox();
                amrex::Array4<Set::Scalar> const &psi = psi_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    psi(i,j,k,0) = 0.0;
                    for (int n = 0; n < material.num_mat; n++)
                        psi(i,j,k,0) += eta(i,j,k,n) * crack.cracktype[n].g_phi(c(i,j,k,0));
                    // psi(i, j, k, 0) = crack.cracktype[0].g_phi(Numeric::Interpolate::NodeToCellAverage(c, i, j, k, 0));
                });
            }
            psi_mf[lev]->FillBoundary();
        }
    }

    void
    TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        // Deciding whether to do an elastic solve or not
        Util::Message(INFO, crack.driving_force_norm, " ", crack.driving_force_reference, " ", crack.driving_force_reference_prev, " ", crack.tol_rel);
        
        // there are basically two types of simulations
        // 1. "Static" simulations where load steps are increased slowly.
        //      In this case, we solve static mechanics problems and quasi-static crack evolution.
        // 2. "Dynamic" simulations where loading rate is specified
        //      In this case, we solve dynamic mechanics problem and quasi-static crack evolution.

        // This decides if the crack solver has converged. If it has, we can do one of the following things
        //      1. Static simulations - Just do elastic solve with no loadstep increase. Clock needs to stay off in this case.
        //      2. Dynamic simulations - Turn on the clock and let the load increase.
        if (    (crack.driving_force_norm / crack.driving_force_reference < crack.tol_rel) ||
                (crack.driving_force_norm < crack.tol_abs) || 
                (crack.crack_prop_iter > crack.max_iter))
        {
            // This means that the crack solver has converged. 
            // Now we need to do elastic solve. 
            crack.crack_prop_iter = 0;
            elastic_do_solve_now = true;

            if (    (m_type == Base::Mechanics<ductile_model>::Type::Static && increase_load_step) || 
                    (m_type == Base::Mechanics<ductile_model>::Type::Dynamic) )
            {
                increase_load_step = false;
                set_new_reference = true;
                RestartClock();
            }
        }

        // This means that the crack field has not converged.
        if (!elastic_do_solve_now)
        {
            crack.crack_prop_iter++;
            StopClock();
            return;
        }

        // Doing an elastic solve
        Util::Message(INFO, "Calling elastic solver");
        Base::Mechanics<ductile_model>::TimeStepBegin(a_time, a_step);

        // Computing pristine energy based on eigen value decomposition
        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());
            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::Array4<Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix> const &stress = stress_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix> const &strain = strain_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const &energy = crack.energy_pristine_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const &history_var = crack.history_var_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    //==========================================================
                    // The formulation below uses decomposition of stress
                    // Set::Matrix sig = stress(i,j,k);
                    // Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(sig);
                    // Set::Vector eValues = eigensolver.eigenvalues();
                    // Set::Matrix eVectors = eigensolver.eigenvectors();

                    // Set::Matrix sig_p = Set::Matrix::Zero();
                    // Set::Matrix sig_n = Set::Matrix::Zero();

                    // for (int n = 0; n < AMREX_SPACEDIM; n++)
                    // {
                    //     if (eValues(n) > 0.0) sig_p += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                    //     else sig_n += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                    // }

                    // Set::Matrix sig_n_dev = sig_n - (1.0/3.0)*sig_n.trace()*Set::Matrix::Identity();

                    // for (int n = 0; n < material.num_mat; n++)
                    // {
                    //     // computing energy based on tensile stress
                    //     energy(i,j,k,0) += eta(i,j,k,n)*material.models[n].W2(sig_p);

                    //     // compute energy from deviatoric portion
                    //     energy(i,j,k,1) += eta(i,j,k,n) * material.models[n].W2(sig_n_dev);
                    // }

                    // // Only update energy if it is increasing. (H^+ according to Miehe)
                    // if (energy(i,j,k,0) < energy_old(i,j,k,0)) energy(i,j,k,0) = energy_old(i,j,k,0);
                    // if (energy(i,j,k,1) < energy_old(i,j,k,1)) energy(i,j,k,1) = energy_old(i,j,k,1);

                    //==========================================================
                    // The formulation below uses Ambati's decomposition of strain
                    // Perform eigenvalue decomposition of strain
                    Set::Matrix eps = strain(i, j, k);
                    Set::Matrix sig = stress(i, j, k);
                    energy(i, j, k, 0) = 0.0;
                    energy(i, j, k, 1) = 0.0;
                    energy(i, j, k, 2) = 0.0;

                    if (!crack.cracktype[0].mixed_mode())
                    {
                        // perform eigenvalue decomposition of strain tensor.
                        Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(eps);
                        Set::Vector eValues = eigensolver.eigenvalues();
                        Set::Matrix eVectors = eigensolver.eigenvectors();

                        // Reconstruct positive and negative counterparts of strain
                        Set::Matrix eps_p = Set::Matrix::Zero();
                        Set::Matrix eps_n = Set::Matrix::Zero();

                        for (int n = 0; n < AMREX_SPACEDIM; n++)
                        {
                            if (eValues(n) > 0.0)
                                eps_p += eValues(n) * (eVectors.col(n) * eVectors.col(n).transpose());
                            else
                                eps_n += eValues(n) * (eVectors.col(n) * eVectors.col(n).transpose());
                        }

                        for (int n = 0; n < material.num_mat; n++){
                            // energy(i, j, k, 0) += eta(i, j, k, n) * material.models[n].W(eps_p);
                            energy(i, j, k, 0) += Numeric::Interpolate::CellToNodeAverage(eta, i, j, k, n) * ( 
                                                    material.models[n].W(eps_p) + material.models[n].PlasticEnergy());
                        }

                        if (energy(i, j, k, 0) > history_var(i, j, k, 0))
                            history_var(i, j, k, 0) = energy(i, j, k, 0);
                    }
                    else
                    {
                        Util::Abort(INFO, "Mixed mode fracture is not supported in ductile fracture");
                    }
                    //=====================================================================
                });
            }
            Util::RealFillBoundary(*crack.energy_pristine_mf[lev], geom[lev]);
            Util::RealFillBoundary(*crack.history_var_mf[lev], geom[lev]);
        }

        integrate_variables_before_advance = false;
        integrate_variables_after_advance = true;
    }

    void
    Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt) override
    {
        if (elastic_do_solve_now)  Base::Mechanics<ductile_model>::Advance(a_lev, a_time, a_dt);

        #if AMREX_SPACEDIM != 2
            Util::Abort(INFO, "This does not work for 1D or 3D yet.");
        #endif

        // advance for crack field
        crack.c_mf[a_lev]->FillBoundary();
        std::swap(crack.c_old_mf[a_lev], crack.c_mf[a_lev]);

        // Create the time integrator
        amrex::TimeIntegrator timeintegrator(*crack.c_old_mf[a_lev], a_time);

        // Set the time integrator RHS - in this case, just relay to our current RHS function
        timeintegrator.set_rhs([&](amrex::MultiFab & rhs_mf, amrex::MultiFab & c_old_mf, const Set::Scalar /*a_time*/)
        {
            // Get the cell size corresponding to this level
            const Set::Scalar* DX = geom[a_lev].CellSize();
            // amrex::Box domain = geom[a_lev].Domain();

            // Iterate over all of the patches on this level
            for (amrex::MFIter mfi(c_old_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                // Get the box (index dimensions) for this patch
                const amrex::Box &bx = mfi.tilebox();

                // Get an array-accessible handle to the data on this patch.
                Set::Patch<const Set::Scalar>  c_old    = c_old_mf.array(mfi);
                Set::Patch<const Set::Scalar>  energy   = crack.history_var_mf[a_lev]->array(mfi);
                Set::Patch<const Set::Scalar>  eta      = material.eta_mf[a_lev]->array(mfi);
                Set::Patch<Set::Scalar>        df       = crack.driving_force_mf[a_lev]->array(mfi);
                Set::Patch<Set::Scalar>        rhs      = rhs_mf.array(mfi);

                // Iterate over the grid on this patch
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {

                    Set::Scalar bilap = 0.0;
                    if (crack.beta > 0.0)
                        bilap = Numeric::Stencil<Set::Scalar, 4, 0, 0>::D(c_old, i, j, k, 0, DX) + Numeric::Stencil<Set::Scalar, 2, 2, 0>::D(c_old, i, j, k, 0, DX) * 2.0 + Numeric::Stencil<Set::Scalar, 0, 4, 0>::D(c_old, i, j, k, 0, DX);
                    
                    Set::Scalar laplacian = Numeric::Laplacian(c_old, i, j, k, 0, DX);
                    if (std::isnan(laplacian)) Util::Message(INFO, "Laplacian is nan at (", i, ", ", j, ")");

                    Set::Scalar Gc = 0.0;
                    Set::Scalar Zeta = 0.0;
                    Set::Scalar Threshold = 0.0;
                    Set::Scalar Mobility = 0.0;
                    Set::Scalar dg_dphi = 0.0;
                    Set::Scalar dw_dphi = 0.0;
                    Set::Scalar _temp_product = 1.0, _temp_product2 = 1.0;

                    for (int m = 0; m < material.num_mat; m++)
                    {
                        Gc          += eta(i, j, k, m) * crack.cracktype[m].Gc(c_old(i, j, k, 0));
                        Zeta        += eta(i, j, k, m) * crack.cracktype[m].Zeta(c_old(i, j, k, 0));
                        Threshold   += eta(i, j, k, m) * crack.cracktype[m].DrivingForceThreshold(c_old(i, j, k, 0));
                        Mobility    += eta(i, j, k, m) * crack.cracktype[m].Mobility(c_old(i, j, k, 0));
                        dg_dphi     += eta(i, j, k, m) * crack.cracktype[m].Dg_phi(c_old(i,j,k));
                        dw_dphi     += eta(i, j, k, m) * crack.cracktype[m].Dw_phi(c_old(i,j,k), 0.0);
                        _temp_product *= eta(i, j, k, m);
                        _temp_product2 *= 0.5;
                    }
                    if (material.num_mat > 1)
                        Gc *= (1.0 - _temp_product * (1. - crack.mult_Gc) / _temp_product2);

                    Set::Scalar en_cell_0 = Numeric::Interpolate::NodeToCellAverage(energy, i, j, k, 0);
                    Set::Scalar en_cell_1 = Numeric::Interpolate::NodeToCellAverage(energy, i, j, k, 1);

                    rhs(i,j,k) = 0.0;


                    if (!crack.cracktype[0].mixed_mode())
                        df(i, j, k, 0) = dg_dphi * (en_cell_0 + en_cell_1) * crack.el_mult / Gc;
                    else
                        df(i, j, k, 0) = dg_dphi * (en_cell_0 + en_cell_1); // * crack.el_mult;
                    rhs(i,j,k) += df(i,j,k,0);

                    df(i, j, k, 1) = dw_dphi / Zeta;
                    rhs(i,j,k) += df(i, j, k, 1);

                    df(i, j, k, 2) = 2.0 * Zeta * laplacian * crack.mult_lap;
                    if (std::isnan(df(i, j, k, 2)))
                        Util::Message(INFO, "DF(2) is nan at (", i, ",", j, ")");
                    rhs(i,j,k) -= df(i, j, k, 2);

                    df(i, j, k, 3) = crack.beta * (0.5 * Zeta * Zeta * Zeta) * bilap;
                    if (std::isnan(df(i, j, k, 3)))
                        Util::Message(INFO, "DF(3) is nan at (", i, ",", j, ")");
                    rhs(i,j,k) += df(i, j, k, 3);

                    df(i, j, k, 4) = std::max(0., rhs(i,j,k) - Threshold);
                    rhs(i,j,k) = - df(i, j, k, 4) * Mobility;
                });
            }
        });

        // We must update the boundaries here and apply boundary conditions 
        timeintegrator.set_post_stage_action([&](amrex::MultiFab & stage_mf, Set::Scalar /*time*/) 
        {
            for (amrex::MFIter mfi(stage_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                Set::Patch<Set::Scalar>  c    = stage_mf.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if (c(i,j,k,0) < 0.0) c(i,j,k,0) = 0.0;
                    else if (c(i,j,k,0) > 1.0) c(i,j,k,0) = 1.0;
                });
            }

            // bc_psi->FillBoundary(stage_mf,0,1,time,0);   
            stage_mf.FillBoundary();

            // bc_df->FillBoundary(*crack.driving_force_mf[a_lev],0,5,time,0);   
            crack.driving_force_mf[a_lev]->FillBoundary();
        });

        // Do the update
        timeintegrator.advance(*crack.c_old_mf[a_lev], *crack.c_mf[a_lev], a_time, a_dt);
    }

    void
    TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        Base::Mechanics<ductile_model>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();

        for (amrex::MFIter mfi(*crack.c_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.convert(amrex::IntVect::TheCellVector());
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(c, i, j, k, 0, DX.data());
                if (grad.lpNorm<2>() * DXnorm > crack.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }

        for (amrex::MFIter mfi(*crack.driving_force_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.convert(amrex::IntVect::TheCellVector());
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &df = crack.driving_force_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(df, i, j, k, 0, DX.data());
                if (grad.lpNorm<2>() * DXnorm > crack.driving_force_refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }

        for (amrex::MFIter mfi(*psi_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.convert(amrex::IntVect::TheCellVector());
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &psi = psi_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(psi, i, j, k, 0, DX.data());
                if (grad.lpNorm<2>() * DXnorm > crack.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }

        if (material.num_mat > 1)
        {
            for (amrex::MFIter mfi(*material.eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.convert(amrex::IntVect::TheCellVector());
                amrex::Array4<char> const &tags = a_tags.array(mfi);
                amrex::Array4<Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix grad = Numeric::Gradient(eta, i, j, k, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > material.m_eta_ref_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                });
            }
        }
    }

    void
    Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &a_box) override
    {
        Set::Vector DX(geom[amrlev].CellSize());
        const Set::Scalar DV = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        amrex::Array4<const Set::Scalar> const &df = crack.driving_force_mf[amrlev]->array(mfi);
        // amrex::Array4<const Set::Scalar> const &c = crack.c_mf[amrlev]->array(mfi); // not used
        // amrex::Array4<const Set::Scalar> const &c_old = crack.c_old_mf[amrlev]->array(mfi); // not used
        amrex::ParallelFor(a_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            crack.driving_force_norm += df(i, j, k, 4) * DV;
        });

        Base::Mechanics<ductile_model>::Integrate(amrlev, time, step, mfi, a_box);
    }

    void
    TimeStepComplete(Set::Scalar /*time*/, int /* iter*/) override
    {
        // if we did an elastic solve, it means that the crack solver converged in TimeStepBegin
        if (elastic_do_solve_now)
            crack.driving_force_reference = crack.driving_force_norm;

        if (m_type == Base::Mechanics<ductile_model>::Type::Static)
        {
            if (crack.driving_force_reference / crack.driving_force_reference_prev < 0.1 || crack.driving_force_reference < 1.e-12)
                increase_load_step = true;

            if (set_new_reference)
            {
                crack.driving_force_reference_prev = crack.driving_force_reference;
                set_new_reference = false;
            }
        }

        if (set_initial_df_once)
        {
            Util::Message(INFO, "First setting the reference value");
            crack.driving_force_reference_prev = crack.driving_force_reference;
            set_initial_df_once = false;
        }

        elastic_do_solve_now = false;
    }

    // Member variables
protected:
    struct
    {
        Set::Field<Set::Scalar> c_mf;
        Set::Field<Set::Scalar> c_old_mf;
        Set::Field<Set::Scalar> energy_pristine_mf;
        Set::Field<Set::Scalar> history_var_mf;
        Set::Field<Set::Scalar> driving_force_mf;
        std::vector<crack_type> cracktype;

        std::vector<std::string> ic_type;
        std::vector<IC::IC<Set::Scalar> *> ic;
        bool is_ic = true;

        Set::Scalar scaleModulusMax = 0.02;
        Set::Scalar refinement_threshold = NAN;
        Set::Scalar mult_lap = NAN;
        Set::Scalar mult_Gc = NAN;
        Set::Scalar el_mult = NAN;

        Set::Scalar driving_force_reference = 1.0;
        Set::Scalar driving_force_reference_prev = 1.0;
        Set::Scalar driving_force_norm = 0.0;
        Set::Scalar driving_force_refinement_threshold = NAN;
        Set::Scalar tol_rel = NAN;
        Set::Scalar tol_abs = NAN;
        Set::Scalar max_iter = NAN;

        Set::Scalar crack_l2_err = 0.0;
        Set::Scalar crack_norm = 0.0;
        Set::Scalar crack_prop_iter = 0;

        Set::Scalar beta = NAN;

        bool plot_df = false;
        bool plot_energy = false;
        bool plot_crack_old = false;
    } crack;
    struct
    {
        Set::Field<Set::Scalar> eta_mf;
        Set::Scalar m_eta_ref_threshold = 0.01;
        std::vector<ductile_model> models;
        IC::IC<Set::Scalar> *ic;
        bool is_ic = false;
        int num_mat = 1;
        bool plot_eta = false;
    } material;

    bool elastic_do_solve_now = true;
    bool increase_load_step = true;
    bool set_initial_df_once = false;
    bool set_new_reference = true;
    int loadstep = 0;

    BC::BC<Set::Scalar> *bc_psi = nullptr;
    BC::BC<Set::Scalar> *bc_df = nullptr;
    BC::Nothing bc_nothing;

    using Base::Mechanics<ductile_model>::m_type;
    using Base::Mechanics<ductile_model>::finest_level;
    using Base::Mechanics<ductile_model>::geom;
    using Base::Mechanics<ductile_model>::model_mf;
    using Base::Mechanics<ductile_model>::psi_mf;
    using Base::Mechanics<ductile_model>::psi_on;
    using Base::Mechanics<ductile_model>::plot_psi;
};
}

#endif
