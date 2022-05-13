///
/// \file PhaseFieldMicrostructure.H
///
#ifndef INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H
#define INTEGRATOR_PHASEFIELDMICROSTRUCTURE_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"
#include "Integrator/Mechanics.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Step.H"
#include "IC/Constant.H"
#include "IC/TabulatedInterface.H"
#include "IC/PerturbedInterface.H"
#include "IC/Voronoi.H"
#include "IC/Sphere.H"
#include "IC/Expression.H"

#include "Model/Interface/GB/GB.H"
#include "Model/Interface/GB/Sin.H"
#include "Model/Interface/GB/AbsSin.H"
#include "Model/Interface/GB/Read.H"
#include "Model/Interface/GB/SH.H"

#include "Model/Solid/Linear/Cubic.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Operator/Elastic.H"

namespace Integrator
{
using model_type = Model::Solid::Affine::Cubic;

enum RegularizationType{
    Wilmore = 1,
    K12 = 2};

///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PhaseFieldMicrostructure : public Integrator // MechanicsBase<Model::Solid::Affine::Cubic>
{
public:
    PhaseFieldMicrostructure() : Integrator() {};
    PhaseFieldMicrostructure(IO::ParmParse &pp) : Integrator() 
    {pp.queryclass(*this);}

    static void Parse(PhaseFieldMicrostructure &value, IO::ParmParse &pp)
    {

        pp.query("pf.number_of_grains", value.number_of_grains);     // Number of grain fields (may be more if using different IC)
        pp.query("pf.M", value.pf.M);                                // Mobility
        if (pp.contains("pf.mu")) pp.query("pf.mu", value.pf.mu);       // Phase field :math:`\mu`
        pp.query("pf.gamma", value.pf.gamma);                        // Phase field :math:`\gamma`
        pp.query("pf.sigma0", value.pf.sigma0);                      // Initial GB energy if not using GB anisotropy
        pp.query("pf.l_gb", value.pf.l_gb);                          // Mobility
        pp.query("pf.elastic_mult",value.pf.elastic_mult);           // Multiplier of elastic energy
        pp.query("pf.elastic_threshold",value.pf.elastic_threshold); // Elastic threshold (:math:`\phi_0`)
        value.pf.L = (4./3.)*value.pf.M / value.pf.l_gb;
    
        pp.query("amr.max_level", value.max_level);            // Maximum AMR level
        pp.query("amr.ref_threshold", value.ref_threshold);    // Phase field refinement threshold

        // Lagrange multiplier method for enforcing volumes
        pp.query("lagrange.on", value.lagrange.on);               // Turn on
        if (value.lagrange.on)
        {
            pp.query("lagrange.lambda", value.lagrange.lambda);   // Lagrange multiplier value
            pp.query("lagrange.vol0", value.lagrange.vol0);       // Prescribed volume
            pp.query("lagrange.tstart", value.lagrange.tstart);   // Time to start enforcing Lagrange multipler
            value.SetThermoInt(1);
        }

        // Anisotropic grain boundary energy parameters
        pp.query("anisotropy.on", value.anisotropy.on);                            // Turn on
        pp.query("anisotropy.beta", value.anisotropy.beta);                        // Regularization param 
        pp.query("anisotropy.tstart", value.anisotropy.tstart);                    // Time to turn on anisotropy
        value.anisotropy.timestep = value.timestep;
        pp.query("anisotropy.timestep", value.anisotropy.timestep);                // Modify timestep when turned on
        value.anisotropy.plot_int = value.plot_int;
        pp.query("anisotropy.plot_int", value.anisotropy.plot_int);                // Modify plot_int when turned on
        value.anisotropy.plot_dt = value.plot_dt;
        pp.query("anisotropy.plot_dt", value.anisotropy.plot_dt);                  // Modify plot_dt when turned on
        pp.query("anisotropy.thermo_int", value.anisotropy.thermo_int);            // Modify thermo int when turned on
        pp.query("anisotropy.thermo_plot_int", value.anisotropy.thermo_plot_int);  // Modify thermo plot int when turned on
        pp.query("anisotropy.elastic_int",value.anisotropy.elastic_int);           // Frequency of elastic calculation

        // Determine the kind of regularization to use
        std::map<std::string, RegularizationType> regularization_type;
        regularization_type["wilmore"] = RegularizationType::Wilmore;
        regularization_type["k12"] = RegularizationType::K12;
        std::string regularization_type_input = "k12";
        pp.query("anisotropy.regularization", regularization_type_input);    // Type of regularization to use  
        value.regularization = regularization_type[regularization_type_input];

        // Set the anisotropic GB model
        pp.query("anisotropy.gb_type", value.gb_type); // Type of GB energy to use
        if      (value.gb_type == "abssin") value.boundary = new Model::Interface::GB::AbsSin(pp,"anisotropy");
        else if (value.gb_type == "sin")    value.boundary = new Model::Interface::GB::Sin(pp,"anisotropy");
        else if (value.gb_type == "read")   value.boundary = new Model::Interface::GB::Read(pp,"anisotropy");
        else if (value.gb_type == "sh")     value.boundary = new Model::Interface::GB::SH(pp,"anisotropy");
        else if (value.anisotropy.on)       Util::Abort(INFO,"A GB model must be specified");


        // Boundary conditions for Etas
        std::string bc_type = "constant";
        pp.query("bc.eta.type",bc_type);  // Type (constnat)
        if (bc_type == "constant")  value.mybc = new BC::Constant(value.number_of_grains,pp,"bc.eta");
        else if (bc_type == "step") value.mybc = new BC::Step(pp,"bc.eta");

        // Initial conditions for Etas
        std::string ic_type = "constant";
        pp.query("ic.type", ic_type); // IC Type
        if (ic_type == "constant")                 value.ic = new IC::Constant(value.geom,pp,"ic");
        else if (ic_type == "perturbed_interface") value.ic = new IC::PerturbedInterface(value.geom,pp,"ic");
        else if (ic_type == "tabulated_interface") value.ic = new IC::TabulatedInterface(value.geom,pp,"ic");
        else if (ic_type == "voronoi")             value.ic = new IC::Voronoi(value.geom, pp, "ic.voronoi");
        else if (ic_type == "expression")          value.ic = new IC::Expression(value.geom,pp,"ic.expression");
        else if (ic_type == "sphere")              value.ic = new IC::Sphere(value.geom,pp); 
        else Util::Abort(INFO, "No valid initial condition specified");


        value.RegisterNewFab(value.eta_new_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "Eta",true);
        value.RegisterNewFab(value.eta_old_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "Eta old",false);

        value.RegisterIntegratedVariable(&value.volume, "volume");
        value.RegisterIntegratedVariable(&value.area, "area");
        value.RegisterIntegratedVariable(&value.gbenergy, "gbenergy");
        value.RegisterIntegratedVariable(&value.realgbenergy, "realgbenergy");
        value.RegisterIntegratedVariable(&value.regenergy, "regenergy");
        value.RegisterIntegratedVariable(&value.elastic.strainenergy, "strainenergy");
        value.RegisterIntegratedVariable(&value.elastic.force, "force");
        value.RegisterIntegratedVariable(&value.elastic.disp, "disp");

        value.RegisterGeneralFab(value.model_mf,1,2);

        // Elasticity
        {
            // Parameters for the elastic solver
            pp.query("elastic.on", value.elastic.on); // Turn on
            if (value.elastic.on)
            {
                value.RegisterNodalFab(value.disp_mf, AMREX_SPACEDIM, 2, "disp",true);
                value.RegisterNodalFab(value.rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
                value.RegisterNodalFab(value.stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
                value.RegisterNodalFab(value.energy_mf, 1, 2, "energy",true);
    
                pp.query("elastic.interval", value.elastic.interval); // How frequently to do elastic solve
                pp.query("elastic.max_coarsening_level", value.elastic.max_coarsening_level); // MLMG max coarsening level
                pp.query("elastic.tol_rel", value.elastic.tol_rel); // Solver relative tolerance
                pp.query("elastic.tol_abs", value.elastic.tol_abs); // Solver absolute tolerance
                pp.query("elastic.tstart", value.elastic.tstart); // Time to start implementing elasticity
    
                pp.queryclass("elastic.bc",value.elastic.bc); // Elastic boundary condition
    
    
                value.elastic.model.resize(value.number_of_grains);
                for (int i = 0; i < value.number_of_grains; i++)
                {
                    model_type mymodel;
                    pp.queryclass("elastic.model",value.elastic.model[i]); // Type of model to use
                }
                pp.queryclass("elastic.model1",value.elastic.model[0]); // Read in for first model
                pp.queryclass("elastic.model2",value.elastic.model[1]); // Read in to second model
            }
        }

    }

    
protected:

    /// \fn    Advance
    /// \brief Evolve phase field in time
    void Advance (int lev, Real time, Real dt) override;
    void Initialize (int lev) override;

    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void TimeStepBegin(amrex::Real time, int iter) override;
    void TimeStepComplete(amrex::Real time, int iter) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &box) override;


    virtual void UpdateModel(int /*a_step*/)
    {
        Util::Abort(INFO,"Not implemented but needs to be.");
    }


private:

    int number_of_grains = 2;
    int number_of_ghost_cells = 3;
    Set::Scalar ref_threshold = 0.1;

    // Cell fab
    Set::Field<Set::Scalar> eta_new_mf; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
    Set::Field<Set::Scalar> eta_old_mf; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
    // Node fab
    Set::Field<Set::Scalar> disp_mf; 
    Set::Field<Set::Scalar> rhs_mf; 
    Set::Field<Set::Scalar> stress_mf;
    Set::Field<Set::Scalar> energy_mf;

    Set::Field<model_type> model_mf;

    BC::BC<Set::Scalar> *mybc;

    //amrex::Real M, mu, gamma, sigma0, l_gb, beta;
    RegularizationType regularization = RegularizationType::K12;
    struct {
        Set::Scalar M;
        Set::Scalar L;
        Set::Scalar mu;
        Set::Scalar gamma;
        Set::Scalar sigma0;
        Set::Scalar l_gb;
        Set::Scalar elastic_mult = 1.0;
        Set::Scalar elastic_threshold = 0.0;
    } pf;

    struct {
        int on = 0;
        Set::Scalar beta;
        Set::Scalar timestep;
        Set::Scalar tstart;
        int plot_int = -1;
        Set::Scalar plot_dt = -1.0;
        int thermo_int = -1, thermo_plot_int = -1;
        Set::Scalar theta0,sigma0,sigma1;
        Set::Scalar phi0 = 0.0;
        int elastic_int = -1;
    } anisotropy;
    
    struct { 
        int on = 0;
        Set::Scalar tstart = 0.0;
        Set::Scalar vol0 = 0.0;
        Set::Scalar lambda = 0.0;
    } lagrange;

    std::string gb_type, filename;

    Model::Interface::GB::GB *boundary;

    IC::IC *ic;

    Set::Scalar volume = 5;
    Set::Scalar area = 0.0;
    Set::Scalar gbenergy = 0.0;
    Set::Scalar realgbenergy = 0.0;
    Set::Scalar regenergy = 0.0;

    struct
    {
        bool        on       = false;
        int         interval = -1;
        int         max_coarsening_level = -1;
        amrex::Real tol_rel = 0.0;
        amrex::Real tol_abs = 1.0E-10;
        amrex::Real tstart = 0.0;
        amrex::Vector<amrex::Real> load_t;
        amrex::Vector<amrex::Real> load_disp;
        std::vector<model_type> model;

        amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xlo,bc_ylo,bc_zlo);
        amrex::Vector<amrex::Real> AMREX_D_DECL(bc_xhi,bc_yhi,bc_zhi);
        std::array<BC::Operator::Elastic::Constant::Type,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xlo, bctype_ylo, bctype_zlo);
        std::array<BC::Operator::Elastic::Constant::Type,AMREX_SPACEDIM> AMREX_D_DECL(bctype_xhi, bctype_yhi, bctype_zhi);
        Operator::Elastic<model_type::sym> op;

        BC::Operator::Elastic::Constant bc;

        Set::Scalar strainenergy = 0.0;
        Set::Scalar force = 0.0;
        Set::Scalar disp = 0.0;
    } elastic;

};
}
#endif
