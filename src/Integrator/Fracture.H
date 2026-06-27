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
#ifndef INTEGRATOR_FRACTURE_H
#define INTEGRATOR_FRACTURE_H

#include "Integrator/Base/Mechanics.H"

#include "Model/Interface/Crack/PFCZM.H"
#include "Model/Solid/Linear/Isotropic.H"

#include "Model/Solid/Finite/NeoHookean.H"

#include "BC/Constant.H"
#include "IC/BMP.H"
#include "IC/Ellipse.H"
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
//using brittle_model = Model::Solid::Linear::Isotropic;
using brittle_model = Model::Solid::Finite::NeoHookean;
using pfczm_crack_type = Model::Interface::Crack::PFCZM;

class Fracture : virtual public Base::Mechanics<brittle_model>
{
public:
    static constexpr const char *name = "fracture";
    Fracture() : Base::Mechanics<brittle_model>() {}

    Fracture(IO::ParmParse &pp) : Base::Mechanics<brittle_model>()
    {
        Parse(*this, pp);
    }

    static void
    Parse(Fracture &value, IO::ParmParse &pp)
    {
        Base::Mechanics<brittle_model>::Parse(value, pp);

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
        pp.queryclass_enumerate<brittle_model>("material.model", value.material.models);
        value.material.num_mat = value.material.models.size();
        Util::Assert(INFO, TEST(value.material.models.size() > 0));

        // This is the fab that stores the material field.
        // For a simple, single material simulation, this will be a uniform field with value 1 everywhere
        value.RegisterNodalFab(value.material.eta_mf, value.material.num_mat, 2, "eta", true);

        // Crack related parsing
        pp.query_default("crack.refinement_threshold", value.crack.refinement_threshold, 0.0001); // mesh refinement criteria
        pp.query_default("crack.df.beta", value.crack.beta, 0.0);                                 // constant multiplier for fourth orther model with bilaplacian
        pp.query_default("crack.df.tol_rel", value.crack.tol_rel, 1E-3);                          // relative tolerance for convergence of driving force
        pp.query_default("crack.df.tol_abs", value.crack.tol_abs, 1E-3);                          // absolute tolerance for convergence of driving force
        pp.query_default("crack.df.max_iter", value.crack.max_iter, 500.);                        // Maximum number of solver iterations
        pp.query_default("crack.df.mult_Gc", value.crack.mult_Gc, 1.0);                           // constant multiplier for controlling fracture energy of interface
        pp.query_default("crack.df.mult_lap", value.crack.mult_lap, 1.0);                         // constant multiplier for second order model with laplacian
        pp.query_default("crack.df.el_mult", value.crack.el_mult, 1.0);                           // unit conversion multiplier between elastic problem and crack problem

        // Let's figure out if there is an initial crack or void
        if (pp.contains("crack.ic.type")) {
            pp.select_enumerate<IC::Constant,IC::Notch,IC::Ellipse, IC::Expression, IC::PNG>("crack.ic", value.crack.ic, value.geom);
            value.crack.is_ic = true;
        }
        else {
            value.crack.is_ic = false;
        }

        value.RegisterNodalFab(value.crack.c_mf, 1, 2, "crack", true);
        value.RegisterNodalFab(value.crack.c_old_mf, 1, 2, "crack_old", false);
        value.RegisterNodalFab(value.crack.driving_force_mf, 5, 2, "driving_force", false);
        value.RegisterNodalFab(value.crack.energy_pristine_mf, 3, 2, "energy_pristine", false);
        value.RegisterNodalFab(value.crack.history_var_mf, 3, 2, "history_variable", false);
        value.RegisterIntegratedVariable(&(value.crack.driving_force_norm), "driving_force_norm");
        value.RegisterIntegratedVariable(&(value.crack.crack_l2_err), "crack_l2_err");
        value.RegisterIntegratedVariable(&(value.crack.crack_norm), "crack_norm");

        // By default the psi variable will be on for fracture simulations
        value.psi_on = true;
        value.bc_psi = new BC::Constant(1, pp, "crack.bc"); // boundary conditions for crack
        value.RegisterNewFab(value.psi_mf, value.bc_psi, 1, 2, "psi", false);

        // This is needed for the specific crack model we are implementing
        pp.queryclass_enumerate<pfczm_crack_type>("crack.model",value.crack.cracktype);
        Util::Assert(INFO, TEST(value.material.models.size() > 0));

        // Needed to store Gc values in a field
        value.RegisterNodalFab(value.visible_Gc_mf,1,2,"Gc",true);
    }

    void
    Initialize(int lev) override
    {
        Base::Mechanics<brittle_model>::Initialize(lev);

        crack.c_mf[lev]->setVal(1.0);
        crack.c_old_mf[lev]->setVal(1.0);

        visible_Gc_mf[lev]->setVal(0.0);

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
        if (m_type == Base::Mechanics<brittle_model>::Type::Disable)
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
                    amrex::Array4<brittle_model> const &model = model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                        model(i, j, k) = brittle_model::Zero();
                        for (int n = 0; n < material.num_mat; n++)
                            model(i, j, k) += eta(i, j, k, n) * material.models[n];
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
                amrex::Array4<const Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    // psi(i,j,k,0) = crack.cracktype[0].g_phi(c(i,j,k,0));
                    psi(i, j, k, 0) = crack.cracktype[0].g_phi(Numeric::Interpolate::NodeToCellAverage(c, i, j, k, 0));
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
        if ((crack.driving_force_norm / crack.driving_force_reference < crack.tol_rel) || crack.crack_prop_iter > crack.max_iter)
        {
            crack.crack_prop_iter = 0;
            elastic_do_solve_now = true;
            increase_load_step = true;
        }
        if ((crack.driving_force_norm < crack.tol_abs) || crack.crack_prop_iter > crack.max_iter)
        {
            elastic_do_solve_now = true;
            crack.crack_prop_iter = 0;
        }
        if (increase_load_step)
        {
            Util::Message(INFO, "Load step = ", loadstep + 1);
            elastic_do_solve_now = true;
            crack.crack_prop_iter = 0;
            increase_load_step = false;
            set_new_reference = true;
            // if (loadstep == 0) set_initial_df_once = true;
            loadstep++;
        }

        // This means that the crack field has not converged.
        if (!elastic_do_solve_now)
        {
            crack.crack_prop_iter++;
            return;
        }

        // Doing an elastic solve
        Base::Mechanics<brittle_model>::TimeStepBegin(a_time, a_step);
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

                    // If brittle_model is finite kinematics, we must convert to eps from F
                    if (brittle_model::kinvar == Model::Solid::F)
                    {
                        eps -= Set::Matrix::Identity();
                        eps = 0.5*(eps + eps.transpose());
                    }
                    
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

                        // If brittle_model is finite kinematics, use F instead of eps.
                        if (brittle_model::kinvar == Model::Solid::F)
                            eps_p += Set::Matrix::Identity();
                    
                        for (int n = 0; n < material.num_mat; n++)
                            energy(i, j, k, 0) += eta(i, j, k, n) * material.models[n].W(eps_p);

                        if (energy(i, j, k, 0) > history_var(i, j, k, 0))
                            history_var(i, j, k, 0) = energy(i, j, k, 0);
                    }
                    else
                    {
                        // Here we implement the 2023 model by Wang et al (DOI: 10.1016/j.apm.2022.12.006)
                        Set::Vector crack_grad = Numeric::Gradient(c, i, j, k, 0, DX.data());
                        if (crack_grad.lpNorm<2>() > 1.e-4)
                        {

                            // perform eigenvalue decomposition of stress tensor.
                            Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(sig);
                            Set::Vector eValues = eigensolver.eigenvalues();
                            Set::Matrix eVectors = eigensolver.eigenvectors();

                            // Let's order the evalues
                            Set::Scalar sig1 = 0., sig2 = 0.;
                            Set::Vector vec1 = Set::Vector::Zero(), vec2 = Set::Vector::Zero();
                            if (eValues(0) > eValues(1))
                            {
                                sig1 = eValues(0);
                                sig2 = eValues(1);
                                vec1 = eVectors.col(0);
                                vec2 = eVectors.col(1);
                            }
                            else
                            {
                                sig1 = eValues(1);
                                sig2 = eValues(0);
                                vec1 = eVectors.col(1);
                                vec2 = eVectors.col(0);
                            }

                            Set::Scalar irwing_length = crack.cracktype[0].l_w();
                            Set::Scalar c_alpha = crack.cracktype[0].c_alpha();

                            Set::Scalar mohr_center = 0.5 * (sig1 + sig2);
                            Set::Scalar mohr_radius = 0.5 * std::abs(sig1 - sig2);

                            Set::Scalar f_theta = 0.0;

                            if (crack.cracktype[0].failure_surface() == Model::Interface::Crack::PFCZM::FSType::WANG2023)
                            {
                                // Storing useful quantities for later.
                                Set::Scalar sig_t = crack.cracktype[0].sig_t();
                                Set::Scalar tau_f = crack.cracktype[0].tau_s();
                                Set::Scalar chi = crack.cracktype[0].chi();
                                Set::Scalar beta_bar = crack.cracktype[0].beta_bar();

                                Set::Scalar beta = 1.0;
                            
                                // Now let's compute a candidate theta
                                // Ideally, theta = 0 is always a valid solution and corresponds to sig_nn = sig1, tau_nm = 0
                                // Set::Scalar theta = 0.0; // removed because: set but not used
                                beta = sig1 < 0 ? beta_bar : 1.0;
                                f_theta = (beta * (crack.el_mult * crack.el_mult * sig1 * sig1 / (sig_t * sig_t))) - 1;
                                
                                // Now we check for other feasible solutions. First we check for beta = 1.
                                Set::Scalar sig_nn1 = 0., tau_nm1 = 0.;
                                Set::Scalar theta_candidate1 = 0.0, f_candidate1 = 0.;
                                Set::Scalar temp = std::abs((mohr_center / mohr_radius) * (chi * chi / (1.0 - (chi * chi))));
                                if (temp >= 0.0 && temp <= 1.0)
                                {
                                    theta_candidate1 = 0.5 * std::acos(temp);

                                    // Now check if this actually leas to a positive sig_nn value
                                    sig_nn1 = mohr_center + (mohr_radius * std::cos(2.0 * theta_candidate1));
                                    tau_nm1 = mohr_radius * std::sin(2.0 * theta_candidate1);

                                    if (sig_nn1 >= 0.0)
                                    {
                                        f_candidate1 = (crack.el_mult * crack.el_mult * sig_nn1 * sig_nn1 / (sig_t * sig_t)) + (crack.el_mult * crack.el_mult * tau_nm1 * tau_nm1 / (tau_f * tau_f)) - 1.0;

                                        // We only use this solution if the value of the failure surface is larger.
                                        if (f_candidate1 > f_theta)
                                        {
                                            f_theta = f_candidate1;
                                            beta = 1.0;
                                        }
                                    }
                                }

                                // Now we check for the third feasible solution for beta = beta_bar.
                                Set::Scalar sig_nn2 = 0., tau_nm2 = 0.;
                                Set::Scalar theta_candidate2 = 0.0, f_candidate2 = 0.;
                                temp = std::abs((mohr_center / mohr_radius) * (beta_bar * chi * chi / (1.0 - (beta_bar * chi * chi))));
                                if (temp >= 0.0 && temp <= 1.0)
                                {
                                    theta_candidate2 = 0.5 * std::acos(temp);

                                    // Now check if this actually leas to a negative sig_nn value
                                    sig_nn2 = mohr_center + (mohr_radius * std::cos(2.0 * theta_candidate2));
                                    tau_nm2 = mohr_radius * std::sin(2.0 * theta_candidate2);

                                    if (sig_nn2 <= 0.0)
                                    {
                                        f_candidate2 = (beta_bar * crack.el_mult * crack.el_mult * sig_nn2 * sig_nn2 / (sig_t * sig_t)) + (crack.el_mult * crack.el_mult * tau_nm2 * tau_nm2 / (tau_f * tau_f)) - 1.0;

                                        // We only use this solution if the value of the failure surface is larger.
                                        if (f_candidate2 > f_theta)
                                        {
                                            f_theta = f_candidate2;
                                            beta = beta_bar;
                                        }
                                    }
                                }
                            }
                            else if (crack.cracktype[0].failure_surface() == Model::Interface::Crack::PFCZM::FSType::MC)
                            {
                                Set::Scalar cohesion = crack.cracktype[0].cohesion();
                                Set::Scalar friction = crack.cracktype[0].friction();
                                // In this case, the angle theta is straightforward.
                                // theta = (pi / 4) - (\phi / 2)
                                Set::Scalar sig_nn = 0., tau_nm = 0.;
                                Set::Scalar theta = std::atan(1) - std::atan(friction);

                                sig_nn = crack.el_mult * (mohr_center + (mohr_radius * std::cos(2.0 * theta)));
                                tau_nm = crack.el_mult * (mohr_radius * std::sin(2.0 * theta));

                                f_theta = (tau_nm + (friction * sig_nn)) * (tau_nm + (friction * sig_nn)) / (cohesion * cohesion);
                                f_theta = f_theta - 1.0;

                                // Util::Message(INFO, "sig1 = ", sig1, ", sig2 = ", sig2, ", f_theta = ", f_theta);
                            }

                            energy(i,j,k,0) = (f_theta > 0) ? (c_alpha / (2.0 * irwing_length)) * ((f_theta)) : 0.0;
                            energy(i,j,k,1) = 0.0;

                            if (energy(i, j, k, 0) + energy(i, j, k, 1) > history_var(i, j, k, 0) + history_var(i, j, k, 1))
                            {
                                history_var(i, j, k, 0) = energy(i, j, k, 0);
                                history_var(i, j, k, 1) = energy(i, j, k, 1);
                            }
                        }
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
        // advance for crack field
        // Util::RealFillBoundary(*crack.c_old_mf[a_lev],geom[a_lev]);
        crack.c_mf[a_lev]->FillBoundary();
        // return;
        std::swap(crack.c_old_mf[a_lev], crack.c_mf[a_lev]);

        const Set::Scalar *DX = geom[a_lev].CellSize();
        amrex::Box domain(geom[a_lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);

        // Trying out the predictor corrector approach.
        //======================================================================
        // Predictor step
        //======================================================================
        for (amrex::MFIter mfi(*crack.c_mf[a_lev], true); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            // bx.grow(1);
            // bx = bx & domain;

            amrex::Array4<const Set::Scalar> const &c_old = crack.c_old_mf[a_lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &energy = crack.history_var_mf[a_lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &eta = material.eta_mf[a_lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &visible_Gc = visible_Gc_mf[a_lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &c = crack.c_mf[a_lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &df = crack.driving_force_mf[a_lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if AMREX_SPACEDIM != 2
                Util::Abort(INFO, "This does not work for 1D or 3D yet.");
#endif

                if (i == lo.x && j == lo.y)
                    c(i, j, k, 0) = c(i + 1, j + 1, k, 0);
                else if (i == lo.x && j == hi.y)
                    c(i, j, k, 0) = c(i + 1, j - 1, k, 0);
                else if (i == hi.x && j == lo.y)
                    c(i, j, k, 0) = c(i - 1, j + 1, k, 0);
                else if (i == hi.x && j == hi.y)
                    c(i, j, k, 0) = c(i - 1, j - 1, k, 0);
                else if (i == lo.x)
                    c(i, j, k) = c(i + 1, j, k, 0);
                else if (j == lo.y)
                    c(i, j, k) = c(i, j + 1, k, 0);
                else if (i == hi.x)
                    c(i, j, k) = c(i - 1, j, k, 0);
                else if (j == hi.y)
                    c(i, j, k) = c(i, j - 1, k, 0);

                // Next, we are doing crack evolution
                else
                {
                    Set::Scalar rhs = 0.0;
                    Set::Scalar bilap = 0.0;
                    if (crack.beta > 0.0)
                    {
                        bilap = Numeric::Stencil<Set::Scalar, 4, 0, 0>::D(c_old, i, j, k, 0, DX) + Numeric::Stencil<Set::Scalar, 2, 2, 0>::D(c_old, i, j, k, 0, DX) * 2.0 + Numeric::Stencil<Set::Scalar, 0, 4, 0>::D(c_old, i, j, k, 0, DX);
                    }

                    // if (std::isnan(bilap))  Util::Message(INFO, "Bilaplacian is nan at (",i,", ",j,")");

                    Set::Scalar laplacian = Numeric::Laplacian(c_old, i, j, k, 0, DX);
                    if (std::isnan(laplacian))
                        Util::Message(INFO, "Laplacian is nan at (", i, ", ", j, ")");

                    Set::Scalar Gc = 0.0;
                    Set::Scalar Zeta = 0.0;
                    Set::Scalar Threshold = 0.0;
                    Set::Scalar Mobility = 0.0;
                    Set::Scalar _temp_product = 1.0, _temp_product2 = 1.0;

                    for (int m = 0; m < material.num_mat; m++)
                    {
                        Gc += eta(i, j, k, m) * crack.cracktype[m].Gc(c_old(i, j, k, 0));
                        Zeta += eta(i, j, k, m) * crack.cracktype[m].Zeta(c_old(i, j, k, 0));
                        Threshold += eta(i, j, k, m) * crack.cracktype[m].DrivingForceThreshold(c_old(i, j, k, 0));
                        Mobility += eta(i, j, k, m) * crack.cracktype[m].Mobility(c_old(i, j, k, 0));
                        _temp_product *= eta(i, j, k, m);
                        _temp_product2 *= 0.5;
                    }
                    if (material.num_mat > 1)
                        Gc *= (1.0 - _temp_product * (1. - crack.mult_Gc) / _temp_product2);

                    visible_Gc(i,j,k) = Gc;

                    // =====================================================================
                    // Set::Scalar en_cell = energy(i, j, k, 0); // set but not used
                    if (!crack.cracktype[0].mixed_mode())
                    {
                        df(i, j, k, 0) = crack.cracktype[0].Dg_phi(c_old(i, j, k)) * (energy(i, j, k, 0) + energy(i, j, k, 1)) * crack.el_mult / Gc;
                    }
                    else
                    {
                        df(i, j, k, 0) = crack.cracktype[0].Dg_phi(c_old(i, j, k)) * (energy(i, j, k, 0) + energy(i, j, k, 1)); // * crack.el_mult;
                    }
                    rhs += df(i, j, k, 0);

                    df(i, j, k, 1) = crack.cracktype[0].Dw_phi(c_old(i, j, k, 0), 0.0) / (Zeta);
                    rhs += df(i, j, k, 1);

                    df(i, j, k, 2) = 2.0 * Zeta * laplacian * crack.mult_lap;
                    if (std::isnan(df(i, j, k, 2)))
                        Util::Message(INFO, "DF(2) is nan at (", i, ",", j, ")");
                    rhs -= df(i, j, k, 2);

                    df(i, j, k, 3) = crack.beta * (0.5 * Zeta * Zeta * Zeta) * bilap;
                    if (std::isnan(df(i, j, k, 3)))
                        Util::Message(INFO, "DF(3) is nan at (", i, ",", j, ")");
                    rhs += df(i, j, k, 3);

                    df(i, j, k, 4) = std::max(0., rhs - Threshold);
                    c(i, j, k, 0) = c_old(i, j, k, 0) - a_dt * df(i, j, k, 4) * Mobility; //*(4.*c_old(i,j,k,0) - 4.*c_old(i,j,k,0)*c_old(i,j,k,0));

                    if (c(i, j, k, 0) < 0.0)
                        c(i, j, k, 0) = 0.0;
                    if (c(i, j, k, 0) > 1.0)
                        c(i, j, k, 0) = 1.0;
                }
            });
        }
        crack.c_mf[a_lev]->FillBoundary();
        crack.driving_force_mf[a_lev]->FillBoundary();

        Base::Mechanics<brittle_model>::Advance(a_lev, a_time, a_dt);
    }

    void
    TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        Base::Mechanics<brittle_model>::TagCellsForRefinement(lev, a_tags, a_time, a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();

        if (!crack.is_ic) {
            // No initial crack
            Set::Scalar stress_refinement_threshold = 1e7;
            for (amrex::MFIter mfi(*stress_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.convert(amrex::IntVect::TheCellVector());
                amrex::Array4<char> const &tags = a_tags.array(mfi);
                amrex::Array4<Set::Matrix> const &stress = stress_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Matrix s = stress(i, j, k);
                    Set::Scalar sxx = s(0,0), syy = s(1,1), sxy = s(0,1);
                    Set::Scalar von_mises = std::sqrt(sxx*sxx - sxx*syy + syy*syy + 3*sxy*sxy);
                    if (von_mises > stress_refinement_threshold)
                        tags(i,j,k) = amrex::TagBox::SET;
                });
            }
        }

        else {
            // Initial crack present
            for (amrex::MFIter mfi(*crack.c_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.convert(amrex::IntVect::TheCellVector());
                amrex::Array4<char> const &tags = a_tags.array(mfi);
                amrex::Array4<Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    Set::Vector grad = Numeric::NodeGradientOnCell(c, i, j, k, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > crack.refinement_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                });
            }
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

        Base::Mechanics<brittle_model>::Integrate(amrlev, time, step, mfi, a_box);
    }

    void Regrid(int lev, Set::Scalar time) override
    {
        if (material.is_ic)
            material.ic->Initialize(lev, material.eta_mf, time);
    }


    void
    TimeStepComplete(Set::Scalar /*time*/, int /* iter*/) override
    {
        if (elastic_do_solve_now)
            crack.driving_force_reference = crack.driving_force_norm;

        if (set_initial_df_once)
        {
            Util::Message(INFO, "First setting the reference value");
            crack.driving_force_reference_prev = crack.driving_force_reference;
            set_initial_df_once = false;
        }

        if (crack.driving_force_reference / (crack.driving_force_reference_prev+1E-8) < 0.1)
        {
            increase_load_step = true;
            crack.driving_force_reference_prev = crack.driving_force_reference;
        }
        else
        {
            Util::Message(INFO,"cur=",crack.driving_force_reference," prev=", crack.driving_force_reference_prev, " ratio=", crack.driving_force_reference / (crack.driving_force_reference_prev+1E-8));
        }
        if (set_new_reference)
        {
            crack.driving_force_reference_prev = crack.driving_force_reference;
            set_new_reference = false;
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

        std::vector<pfczm_crack_type> cracktype;

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
    } crack;
    struct
    {
        Set::Field<Set::Scalar> eta_mf;
        Set::Scalar m_eta_ref_threshold = 0.01;
        std::vector<brittle_model> models;
        IC::IC<Set::Scalar> *ic;
        bool is_ic = false;
        int num_mat = 1;
    } material;

    bool elastic_do_solve_now = true;
    bool increase_load_step = true;
    bool set_initial_df_once = false;
    bool set_new_reference = true;
    int loadstep = 0;

    BC::BC<Set::Scalar> *bc_psi = nullptr;

    Set::Field<Set::Scalar>visible_Gc_mf;

    using Base::Mechanics<brittle_model>::m_type;
    using Base::Mechanics<brittle_model>::finest_level;
    using Base::Mechanics<brittle_model>::geom;
    using Base::Mechanics<brittle_model>::model_mf;
    using Base::Mechanics<brittle_model>::psi_mf;
    using Base::Mechanics<brittle_model>::psi_on;
};
}

#endif
