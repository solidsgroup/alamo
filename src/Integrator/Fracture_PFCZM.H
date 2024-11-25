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
#ifndef INTEGRATOR_FRACTURE_PFCZM_H
#define INTEGRATOR_FRACTURE_PFCZM_H

#include "Integrator/Base/Mechanics.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Interface/Crack/PFCZM.H"

#include "IC/IC.H"
#include "IC/Ellipsoid.H"
#include "IC/Notch.H"
#include "IC/Laminate.H"
#include "IC/PerturbedInterface.H"
#include "IC/BMP.H"
#include "IC/Expression.H"

#include "Numeric/Stencil.H"

#include <eigen3/Eigen/Dense>

namespace Integrator
{
using brittle_model = Model::Solid::Linear::Isotropic;
using pfczm_crack_type = Model::Interface::Crack::PFCZM;

class Fracture_PFCZM: 
    virtual public Base::Mechanics<brittle_model>
{

public:
    Fracture_PFCZM(): Base::Mechanics<brittle_model>(){}

    Fracture_PFCZM(IO::ParmParse &pp) : Base::Mechanics<brittle_model>()
    {
        Parse(*this,pp);
    }

    static void Parse(Fracture_PFCZM &value, IO::ParmParse &pp)
    {
        Base::Mechanics<brittle_model>::Parse(value,pp);

        // Material field related parsing
        pp.query("material.refinement_threshold", value.material.m_eta_ref_threshold);

        // Let's figure out if we are doing multi-material simulations.
        std::string mat_ic_type;
        pp.query("material.ic.type", mat_ic_type); // type of material - single material vs multimaterial
        if(mat_ic_type == "laminate")
        {
            IC::Laminate *tmpic = new IC::Laminate(value.geom);
            pp.queryclass("material.ic.laminate", *tmpic);  // get parameters of a laminate inclusion
            pp.query("material.ic.laminate.number_of_inclusions", value.material.num_mat); // number of laminate inclusions
            value.material.num_mat++;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else if(mat_ic_type == "ellipse")
        {
            IC::Ellipse *tmpic = new IC::Ellipse(value.geom);
            pp.queryclass("material.ic.ellipse", *tmpic);   // get parameters of elliptical inclusion
            value.material.num_mat = 2;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else if(mat_ic_type == "perturbed_interface")
        {
            IC::PerturbedInterface *tmpic = new IC::PerturbedInterface(value.geom);
            pp.queryclass("material.ic.perturbed_interface", *tmpic);   // get parameters of a wavy interface composite
            value.material.num_mat = 2;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else if(mat_ic_type == "bmp")
        {
            IC::BMP *tmpic = new IC::BMP(value.geom);
            pp.queryclass("material.ic.bmp", *tmpic);   // get parameters for image based material initialization
            value.material.num_mat = 2;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        // This is the simplest case - homogeneous material
        else
        {
            value.material.num_mat = 1;
            value.material.is_ic = false;
        }

        // Let's query material properties for different materials in the system
        for (int i = 0; i<value.material.num_mat; i++)
        {
            std::string name = "material.model" + std::to_string(i+1);
            std::string name_type = "material.model" + std::to_string(i+1) + ".type";
            std::string mat_model_type;
            pp.query(name_type.c_str(),mat_model_type); // read type of material - isotropic is the only option for now.
            
            // For now we only have isotropic implemented. 
            if(mat_model_type == "isotropic")
            {
                brittle_model tmp_model;
                pp.queryclass(std::string((name+".isotropic").data()), tmp_model);  // read values of isotropic model
                value.material.models.push_back(tmp_model);
            }
            else
                Util::Abort(INFO, "This type of material model is not supported for fracture.");
        }
        Util::Assert(INFO,TEST(value.material.models.size() > 0));

        // This is the fab that stores the material field. 
        // For a simple, single material simulation, this will be a uniform field with value 1 everywhere
        value.RegisterNodalFab(value.material.eta_mf, value.material.num_mat, 2, "eta", true);

        // Crack related parsing
        pp.query("crack.refinement_threshold",value.crack.refinement_threshold); // mesh refinement criteria
        pp.query("crack.beta", value.crack.beta);   // constant multiplier for fourth orther model with bilaplacian
        pp.query("crack.df.tol_rel",value.crack.tol_rel); // relative tolerance for convergence of driving force
        pp.query("crack.df.tol_abs",value.crack.tol_abs);   // absolute tolerance for convergence of driving force
        pp.query("crack.df.mult_Gc", value.crack.mult_Gc);  // constant multiplier for controlling fracture energy of interface
        pp.query("crack.df.mult_lap", value.crack.mult_lap);    // constant multiplier for second order model with laplacian
        pp.query("crack.df.el_mult", value.crack.el_mult);  // unit conversion multiplier between elastic problem and crack problem

        // Let's figure out if there is an initial crack or void
        pp.queryarr("crack.ic.type",value.crack.ic_type);   // IC for initial crack
        if(value.crack.ic_type.size() == 0) value.crack.is_ic = false;
        else
        {
            value.crack.is_ic= true;
            for (unsigned int i = 0; i<value.crack.ic_type.size(); i++)
            {
                if(value.crack.ic_type[i] == "notch")
                {
                    IC::Notch *tmpic = new IC::Notch(value.geom);
                    pp.queryclass("crack.ic.notch",*tmpic); // get parameters of notch type crack
                    value.crack.ic.push_back(tmpic);
                }
                else if(value.crack.ic_type[i] == "ellipsoid")
                {
                    IC::Ellipsoid *tmpic = new IC::Ellipsoid(value.geom);
                    pp.queryclass("crack.ic.ellipsoid",*tmpic); // get parameters of ellipsoidal crack
                    value.crack.ic.push_back(tmpic);
                }
                else if(value.crack.ic_type[i] == "expression")
                {
                    IC::Expression *tmpic = new IC::Expression(value.geom);
                    pp.queryclass("crack.ic.expression",*tmpic); // get parameters of the expression crack
                    value.crack.ic.push_back(tmpic);
                }
                else
                    Util::Abort(INFO, "This IC hasn't been implemented yet");
            }
        }
        
        value.RegisterNodalFab(value.crack.c_mf, 1, 2, "crack", true);
        value.RegisterNodalFab(value.crack.c_old_mf, 1, 2, "crack_old", true);
        value.RegisterNodalFab(value.crack.driving_force_mf, 5, 2, "driving_force", true);
        value.RegisterNodalFab(value.crack.energy_pristine_mf, 3, 2, "energy_pristine", true);
        value.RegisterNodalFab(value.crack.energy_pristine_old_mf, 3, 2, "energy_pristine_old", false);
        value.RegisterIntegratedVariable(&(value.crack.driving_force_norm),"driving_force_norm");
        value.RegisterIntegratedVariable(&(value.crack.crack_l2_err), "crack_l2_err");
        value.RegisterIntegratedVariable(&(value.crack.crack_norm), "crack_norm");
        
        // By default the psi variable will be on for fracture simulations
        value.psi_on = true;
        value.bc_psi = new BC::Constant(1,pp,"crack.bc");   // boundary conditions for crack
        value.RegisterNewFab(value.psi_mf,value.bc_psi,1,2,"psi",true);
        
        // This is needed for the specific crack model we are implementing
        for (int i = 0; i<value.material.num_mat; i++)
        {
            std::string name = "crack.model" + std::to_string(i+1);
            std::string name_type = "crack.model" + std::to_string(i+1) + ".type";
            std::string mat_model_type;
            pp.query(name_type.c_str(),mat_model_type); // type of crack model - only constant so far
            // For now we only have a crack of constant Gc implemented
            if(mat_model_type == "pfczm")
            {
                pfczm_crack_type tmp_model;
                pp.queryclass(std::string((name+".pfczm").data()), tmp_model);   // get parameters for a constant crack model.
                value.crack.cracktype.push_back(tmp_model);
            }
            else
                Util::Abort(INFO, "This type of material model is not supported for fracture.");
        }
        Util::Assert(INFO,TEST(value.material.models.size() > 0));
    }

    void Initialize(int lev) override
    {
        Base::Mechanics<brittle_model>::Initialize(lev);

        crack.c_mf[lev]->setVal(1.0);
        crack.c_old_mf[lev]->setVal(1.0);
        
        // For now, we are setting psi field to 1. This prevents nans.
        // We will modify the psi field later.
        psi_mf[lev]->setVal(1.0);

        if(crack.is_ic)
        {
            for (int i=0; i < crack.ic.size(); i++)
            {
                // This initializes the crack field.
                crack.ic[i]->Add(lev,crack.c_mf);
                crack.ic[i]->Add(lev,crack.c_old_mf);
            }
        }

        // This initializes the material field
        material.eta_mf[lev]->setVal(1.0);
        if (material.is_ic) material.ic->Initialize(lev,material.eta_mf);
        else material.eta_mf[lev]->setVal(1.0);
    }

    virtual void UpdateModel(int a_step, Set::Scalar /*a_time*/) override
    {
        if (m_type == Base::Mechanics<brittle_model>::Type::Disable) return;

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
                                                model(i,j,k) = brittle_model::Zero();
                                                for (int n = 0; n < material.num_mat; n++) 
                                                    model(i,j,k) += eta(i,j,k,n) * material.models[n];
                                            });
                }
                Util::RealFillBoundary(*model_mf[lev],geom[lev]);
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

                                            psi(i,j,k,0) = crack.cracktype[0].g_phi(Numeric::Interpolate::NodeToCellAverage(c,i,j,k,0));
                                        });
            }
            psi_mf[lev]->FillBoundary();
        }
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        // Deciding whether to do an elastic solve or not
        Util::Message(INFO,crack.driving_force_norm," ",crack.driving_force_reference," ",crack.tol_rel);
        if ((crack.driving_force_norm / crack.driving_force_reference < crack.tol_rel) || crack.crack_prop_iter > 100)
        {
            crack.crack_prop_iter = 0;
            elastic_do_solve_now = true;
        }
        if ((crack.driving_force_norm < crack.tol_abs) || crack.crack_prop_iter > 100)
        {
           elastic_do_solve_now = true;
           crack.crack_prop_iter = 0;
        }
        
        // This means that the crack field has not converged.
        if (!elastic_do_solve_now) {crack.crack_prop_iter++; return;}

        // Doing an elastic solve
        Base::Mechanics<brittle_model>::TimeStepBegin(a_time,a_step);

        // Computing pristine energy based on eigen value decomposition
        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            std::swap(crack.energy_pristine_mf[lev],crack.energy_pristine_old_mf[lev]);
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());

            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox(); // & domain;
                amrex::Array4<Set::Scalar>     const &eta   = material.eta_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>     const &strain  = strain_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>     const &energy  = crack.energy_pristine_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>     const &energy_old  = crack.energy_pristine_old_mf[lev]->array(mfi);
                
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                {
                    // Perform eigenvalue decomposition of strain
                    Set::Matrix eps = strain(i,j,k);

                    // We can also do a volumetric deviatoric decomposition.
                    Set::Matrix eps_sph = (1.0/3.0) * eps.trace() * Set::Matrix::Identity();
                    Set::Matrix eps_dev = eps - eps_sph;

                    Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(eps);
                    Set::Vector eValues = eigensolver.eigenvalues();
                    Set::Matrix eVectors = eigensolver.eigenvectors();

                    // Reconstruct positive and negative counterparts of strain
                    Set::Matrix eps_p = Set::Matrix::Zero();
                    Set::Matrix eps_n = Set::Matrix::Zero();

                    for (int n = 0; n<AMREX_SPACEDIM; n++)
                    {
                        if(eValues(n) > 0.0) eps_p += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                        else eps_n += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                    }

                    for (int n = 0; n < material.num_mat; n++)
                    {
                        // computing energy based on Miehe
                        energy(i,j,k,0) += eta(i,j,k,n)*material.models[n].W(eps_p);

                        // computing energy from volumetric portion.
                        if (eps.trace() >= 0) energy(i,j,k,1) += eta(i,j,k,n)*material.models[n].W(eps_sph);
                        
                        // compute energy from deviatoric portion
                        energy(i,j,k,2) += eta(i,j,k,n) * material.models[n].W(eps_dev);
                    }
                    
                    // Only update energy if it is increasing. (H^+ according to Miehe)
                    if (energy(i,j,k,0) < energy_old(i,j,k,0)) energy(i,j,k,0) = energy_old(i,j,k,0);
                    if (energy(i,j,k,1) < energy_old(i,j,k,1)) energy(i,j,k,1) = energy_old(i,j,k,1);
                    if (energy(i,j,k,2) < energy_old(i,j,k,2)) energy(i,j,k,2) = energy_old(i,j,k,2);
                });
            }
            Util::RealFillBoundary(*crack.energy_pristine_mf[lev],geom[lev]);
        }

        integrate_variables_before_advance = false;
        integrate_variables_after_advance = true;
    }

    void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt)
    {
        // advance for crack field
        // Util::RealFillBoundary(*crack.c_old_mf[a_lev],geom[a_lev]);
        crack.c_mf[a_lev]->FillBoundary();
        // return;
        std::swap (crack.c_old_mf[a_lev], crack.c_mf[a_lev]);

        const Set::Scalar* DX = geom[a_lev].CellSize();
        amrex::Box domain(geom[a_lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

        // Trying out the predictor corrector approach.
        //======================================================================
        // Predictor step
        //======================================================================
        for (amrex::MFIter mfi(*crack.c_mf[a_lev], true); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            // bx.grow(1);
            // bx = bx & domain;

            amrex::Array4<const Set::Scalar> const& c_old = crack.c_old_mf[a_lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& energy = crack.energy_pristine_mf[a_lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& eta = material.eta_mf[a_lev]->array(mfi);

            amrex::Array4<Set::Scalar> const& c = crack.c_mf[a_lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& df = crack.driving_force_mf[a_lev]->array(mfi);

            amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
                #if AMREX_SPACEDIM !=2
                Util::Abort(INFO, "This does not work for 1D or 3D yet.");
                #endif

                if      (i == lo.x && j == lo.y) c(i,j,k,0) = c(i+1,j+1,k,0);
                else if (i == lo.x && j == hi.y) c(i,j,k,0) = c(i+1,j-1,k,0);
                else if (i == hi.x && j == lo.y) c(i,j,k,0) = c(i-1,j+1,k,0);
                else if (i == hi.x && j == hi.y) c(i,j,k,0) = c(i-1,j-1,k,0);
                else if (i == lo.x) c(i,j,k) = c(i+1,j,k,0);
                else if (j == lo.y) c(i,j,k) = c(i,j+1,k,0);
                else if (i == hi.x) c(i,j,k) = c(i-1,j,k,0);
                else if (j == hi.y) c(i,j,k) = c(i,j-1,k,0);

                // Next, we are doing crack evolution
                else
                {
                    Set::Scalar rhs = 0.0;
                    // Set::Scalar bilap = 0.0;
                    // if (crack.beta > 0.0)
                    // {
                    //     bilap = 
                    //         Numeric::Stencil<Set::Scalar, 4, 0, 0>::D(c_old,i,j,k,0,DX) +
                    //         Numeric::Stencil<Set::Scalar, 2, 2, 0>::D(c_old,i,j,k,0,DX) * 2.0 +
                    //         Numeric::Stencil<Set::Scalar, 0, 4, 0>::D(c_old,i,j,k,0,DX);
                    // }

                    // if (std::isnan(bilap))  Util::Message(INFO, "Bilaplacian is nan at (",i,", ",j,")");
                    
                    Set::Scalar laplacian = Numeric::Laplacian(c_old,i,j,k,0,DX);
                    if (std::isnan(laplacian))  Util::Message(INFO, "Laplacian is nan at (",i,", ",j,")");

                    Set::Scalar Gc = 0.0;
                    Set::Scalar Zeta = 0.0;
                    Set::Scalar Threshold = 0.0;
                    Set::Scalar Mobility = 0.0;
                    Set::Scalar _temp_product = 1.0, _temp_product2 = 1.0;
                    
                    for (int m = 0; m < material.num_mat; m++)
                    {
                        Gc += eta(i,j,k,m)*crack.cracktype[m].Gc(c_old(i,j,k,0));
                        Zeta += eta(i,j,k,m)*crack.cracktype[m].Zeta(c_old(i,j,k,0));
                        Threshold += eta(i,j,k,m)*crack.cracktype[m].DrivingForceThreshold(c_old(i,j,k,0));
                        Mobility += eta(i,j,k,m)*crack.cracktype[m].Mobility(c_old(i,j,k,0));
                        _temp_product *= eta(i,j,k,m);
                        _temp_product2 *= 0.5;
                    }
                    if (material.num_mat > 1) Gc *= (1.0 - _temp_product*(1.-crack.mult_Gc)/_temp_product2);

                    // =====================================================================
                    Set::Scalar en_cell = energy(i,j,k,0);
                    if (! crack.cracktype[0].mixed_mode())
                        df(i,j,k,0) = crack.cracktype[0].Dg2_phi(c_old(i,j,k))*en_cell*crack.el_mult / Gc;
                    else 
                    {
                        Set::Scalar GcI = crack.cracktype[0].GcI(0);
                        Set::Scalar GcII = crack.cracktype[0].GcII(0);
                        df(i,j,k,0) = crack.cracktype[0].Dg2_phi(c_old(i,j,k))*crack.el_mult*((energy(i,j,k,1)/GcI) + (energy(i,j,k,2)/GcII));
                        df(i,j,k,0) *= crack.cracktype[0].correction_eta();
                    }
                    rhs += df(i,j,k,0);

                    df(i,j,k,1) = crack.cracktype[0].Dw_phi(c_old(i,j,k,0),0.0)/(Zeta);
                    rhs += df(i,j,k,1);

                    df(i,j,k,2) = 2.0*Zeta*laplacian * crack.mult_lap;
                    if (std::isnan(df(i,j,k,2))) Util::Message(INFO, "DF(2) is nan at (",i,",",j,")");
                    rhs -= df(i,j,k,2);

                    // df(i,j,k,3) = crack.beta*(0.5*Zeta*Zeta*Zeta)*bilap;
                    // if (std::isnan(df(i,j,k,3))) Util::Message(INFO, "DF(3) is nan at (",i,",",j,")");
                    // rhs += df(i,j,k,3);
                   
                    df(i,j,k,4) = std::max(0.,rhs - Threshold);
                    c(i,j,k,0) = c_old(i,j,k,0) - a_dt*df(i,j,k,4)*Mobility; //*(4.*c_old(i,j,k,0) - 4.*c_old(i,j,k,0)*c_old(i,j,k,0));

                    if (c (i,j,k,0) < 0.0) c(i,j,k,0) = 0.0;
                    if (c (i,j,k,0) > 1.0) c(i,j,k,0) = 1.0;
                }
            });
        }
        crack.c_mf[a_lev]->FillBoundary();
        crack.driving_force_mf[a_lev]->FillBoundary();

        Base::Mechanics<brittle_model>::Advance(a_lev,a_time,a_dt);
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        Base::Mechanics<brittle_model>::TagCellsForRefinement(lev,a_tags,a_time,a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();

        for (amrex::MFIter mfi(*crack.c_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.convert(amrex::IntVect::TheCellVector());
            amrex::Array4<char> const& tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector grad = Numeric::NodeGradientOnCell(c, i, j, k, DX.data());
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
                amrex::Array4<char> const& tags = a_tags.array(mfi);
                amrex::Array4<Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Matrix grad = Numeric::Gradient(eta, i, j, k, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > material.m_eta_ref_threshold)
                            tags(i, j, k) = amrex::TagBox::SET;
                });
            }
        }
    }

    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &a_box) override
    {
        Set::Vector DX(geom[amrlev].CellSize());
        const Set::Scalar DV = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        amrex::Array4<const Set::Scalar> const &df = crack.driving_force_mf[amrlev]->array(mfi);
        amrex::Array4<const Set::Scalar> const &c = crack.c_mf[amrlev]->array(mfi);
        amrex::Array4<const Set::Scalar> const &c_old = crack.c_old_mf[amrlev]->array(mfi);
        amrex::ParallelFor(a_box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            crack.driving_force_norm += df(i,j,k,4)*DV;
        });

        Base::Mechanics<brittle_model>::Integrate(amrlev,time,step,mfi,a_box);
    }

    void TimeStepComplete(Set::Scalar /*time*/, int /* iter*/)
    {
        if (elastic_do_solve_now)
            crack.driving_force_reference = crack.driving_force_norm;
        elastic_do_solve_now = false;
    }
    

// Member variables
protected:
    struct{
        Set::Field<Set::Scalar> c_mf;
        Set::Field<Set::Scalar> c_old_mf;
        // Set::Field<Set::Scalar> c_mid_mf;
        Set::Field<Set::Scalar> energy_pristine_mf;
        Set::Field<Set::Scalar> energy_pristine_old_mf;
        Set::Field<Set::Scalar> driving_force_mf;
        std::vector<pfczm_crack_type> cracktype;

        std::vector<std::string> ic_type;
        std::vector<IC::IC*> ic;
        bool is_ic = true;

        Set::Scalar scaleModulusMax = 0.02;
        Set::Scalar refinement_threshold = 0.0001;
        Set::Scalar mult_lap = 1.0;
        Set::Scalar mult_Gc = 1.0;
        Set::Scalar el_mult = 1.0;

        Set::Scalar driving_force_reference = 1.0;
        Set::Scalar driving_force_norm = 0.0;
        Set::Scalar tol_rel = 1.e-3;
        Set::Scalar tol_abs = 1.e-3;
        
        Set::Scalar crack_l2_err = 0.0;
        Set::Scalar crack_norm = 0.0;
        Set::Scalar crack_prop_iter = 0;

        Set::Scalar beta = 0.0;
    } crack;
    struct{
        Set::Field<Set::Scalar> eta_mf;
        Set::Scalar m_eta_ref_threshold = 0.01;
        std::vector<brittle_model> models;
        IC::IC *ic;
        bool is_ic = false;
        int num_mat = 1;
    } material;

    bool elastic_do_solve_now = true;
    BC::BC<Set::Scalar> *bc_psi = nullptr;

    using Base::Mechanics<brittle_model>::m_type;
    using Base::Mechanics<brittle_model>::finest_level;
    using Base::Mechanics<brittle_model>::geom;
    using Base::Mechanics<brittle_model>::model_mf;
    using Base::Mechanics<brittle_model>::psi_mf;
    using Base::Mechanics<brittle_model>::psi_on;
};
}


#endif
