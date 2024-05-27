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
#include "IC/Ellipse.H"

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
//using model_type = Model::Solid::Affine::Cubic;

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
template<class model_type>
class PhaseFieldMicrostructure : public Base::Mechanics<model_type>
{
public:
    PhaseFieldMicrostructure() : Integrator() {};
    PhaseFieldMicrostructure(IO::ParmParse &pp) : Integrator() 
    {pp.queryclass(*this);}

    static void Parse(PhaseFieldMicrostructure &value, IO::ParmParse &pp)
    {
        BL_PROFILE("PhaseFieldMicrostructure::Parse");
        
        pp.query("pf.number_of_grains", value.number_of_grains);     // Number of grain fields (may be more if using different IC)
        pp.query("pf.M", value.pf.M);                                // Mobility
        if (pp.contains("pf.mu")) pp.query("pf.mu", value.pf.mu);       // Phase field :math:`\mu`
        pp.query("pf.gamma", value.pf.gamma);                        // Phase field :math:`\gamma`
        pp.query("pf.sigma0", value.pf.sigma0);                      // Initial GB energy if not using GB anisotropy
        pp.query("pf.l_gb", value.pf.l_gb);                          // Mobility
        pp.query("pf.elastic_df",value.pf.elastic_df);               // Determine whether to use elastic driving force
        pp.query("pf.elastic_mult",value.pf.elastic_mult);           // Multiplier of elastic energy

        pp.query("pf.threshold.value",value.pf.threshold.value);        // Value used for thresholding kinetic relation
        pp.query("pf.threshold.chempot",value.pf.threshold.chempot);    // Whether to include chemical potential in threshold
        pp.query("pf.threshold.boundary",value.pf.threshold.boundary);  // Whether to include boundary energy in threshold
        pp.query("pf.threshold.corner",value.pf.threshold.corner);      // Whether to include corner regularization in threshold
        pp.query("pf.threshold.lagrange",value.pf.threshold.lagrange);  // Whether to include lagrange multiplier in threshold
        pp.query("pf.threshold.mechanics",value.pf.threshold.mechanics);// Whether to include mechanical driving force in threshold
        {
            std::string type = "continuous";
            pp.query("pf.threshold.type",type); // Type of thresholding to use
            if (type == "continuous") value.pf.threshold.type = ThresholdType::Continuous;
            else if (type == "chop") value.pf.threshold.type = ThresholdType::Chop;
            else Util::Abort(INFO,"Invalid threshold type");
        }
        
        value.pf.L = (4./3.)*value.pf.M / value.pf.l_gb;
    
        pp.query("amr.max_level", value.max_level);            // Maximum AMR level
        pp.query("amr.ref_threshold", value.ref_threshold);    // Phase field refinement threshold

        // Elasticity
        pp.query("mechanics.tstart",value.mechanics.tstart);
        value.mechanics.model.clear();
        value.mechanics.model.resize(value.number_of_grains,model_type::Zero());
        for (int i=0; i < value.number_of_grains; i++)
        {
            // By default, read in the model specified by "mechanics.model"
            if (pp.getEntries("mechanics.model").size())
                pp.queryclass("mechanics.model", value.mechanics.model[i]);

            // If individual models are specified, read those in and overwrite
            std::string name = "mechanics.model" + std::to_string(i+1);
            if (pp.getEntries(name).size())
                pp.queryclass(std::string(name.data()), value.mechanics.model[i]);
        }

        value.m_type = Base::Mechanics<model_type>::Disable; // Turn mechanics off by default
        pp.queryclass<Base::Mechanics<model_type>>("mechanics",&value);
        if (value.m_type == Base::Mechanics<model_type>::Type::Static) 
            value.number_of_ghost_cells = std::max(value.number_of_ghost_cells,3);


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
        if (value.anisotropy.on) 
            value.number_of_ghost_cells = std::max(value.number_of_ghost_cells,2);

        // Determine the kind of regularization to use
        std::map<std::string, RegularizationType> regularization_type;
        regularization_type["wilmore"] = RegularizationType::Wilmore;
        regularization_type["k12"] = RegularizationType::K12;
        std::string regularization_type_input = "k12";
        pp.query("anisotropy.regularization", regularization_type_input);    // Type of regularization to use  
        value.regularization = regularization_type[regularization_type_input];

        // Set the anisotropic GB model
        pp.query("anisotropy.gb_type", value.gb_type); // Type of GB energy to use
        if      (value.gb_type == "abssin") value.boundary = new Model::Interface::GB::AbsSin(pp,"anisotropy.abssin");
        else if (value.gb_type == "sin")    value.boundary = new Model::Interface::GB::Sin(pp,"anisotropy.sin");
        else if (value.gb_type == "read")   value.boundary = new Model::Interface::GB::Read(pp,"anisotropy.read");
        else if (value.gb_type == "sh")     value.boundary = new Model::Interface::GB::SH(pp,"anisotropy.sh");
        else if (value.anisotropy.on)       Util::Abort(INFO,"Invalid gb type " + value.gb_type);

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
        else if (ic_type == "sphere")              value.ic = new IC::Sphere(value.geom,pp,"ic.sphere"); 
        else if (ic_type == "ellipse")             value.ic = new IC::Ellipse(value.geom,pp,"ic.ellipse");
        else Util::Abort(INFO, "No valid initial condition specified");

        // Anisotropic mobility
        pp.query("anisotropic_kinetics.on",value.anisotropic_kinetics.on);
        if (value.anisotropic_kinetics.on)
        {
            pp.query("anisotropic_kinetics.tstart",value.anisotropic_kinetics.tstart);
            std::string mobility_filename, threshold_filename;
            pp.query("anisotropic_kinetics.mobility",mobility_filename);
            value.anisotropic_kinetics.mobility = Numeric::Interpolator::Linear<Set::Scalar>::Read(mobility_filename);
            pp.query("anisotropic_kinetics.threshold",threshold_filename);
            value.anisotropic_kinetics.threshold = Numeric::Interpolator::Linear<Set::Scalar>::Read(threshold_filename);
            value.RegisterNewFab(value.anisotropic_kinetics.L_mf, value.mybc, value.number_of_grains, 0, "mobility",true);
            value.RegisterNewFab(value.anisotropic_kinetics.threshold_mf, value.mybc, value.number_of_grains, 0, "theshold",true);
        }




        value.RegisterNewFab(value.eta_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "Eta",true);
        value.RegisterNewFab(value.driving_force_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "DrivingForce",true);
        value.RegisterNewFab(value.driving_force_threshold_mf, value.mybc, value.number_of_grains, value.number_of_ghost_cells, "DrivingForceThreshold",true);

        value.RegisterIntegratedVariable(&value.volume, "volume");
        value.RegisterIntegratedVariable(&value.area, "area");
        value.RegisterIntegratedVariable(&value.gbenergy, "gbenergy");
        value.RegisterIntegratedVariable(&value.realgbenergy, "realgbenergy");
        value.RegisterIntegratedVariable(&value.regenergy, "regenergy");


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


    virtual void UpdateModel(int /*a_step*/) override;

private:

    int number_of_grains = 2;
    int number_of_ghost_cells = 1;
    Set::Scalar ref_threshold = 0.1;

    // Cell fab
    Set::Field<Set::Scalar> eta_mf; // Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
    Set::Field<Set::Scalar> driving_force_mf;
    Set::Field<Set::Scalar> driving_force_threshold_mf;

    BC::BC<Set::Scalar> *mybc;

    //amrex::Real M, mu, gamma, sigma0, l_gb, beta;
    RegularizationType regularization = RegularizationType::K12;
    enum ThresholdType {
        Continuous = 0, Chop = 1
    };
    struct {
        Set::Scalar M;
        Set::Scalar L;
        Set::Scalar mu;
        Set::Scalar gamma;
        Set::Scalar sigma0;
        Set::Scalar l_gb;
        bool elastic_df = false;
        Set::Scalar elastic_mult = 1.0;
        struct {
            bool chempot = true;
            bool boundary = true;
            bool corner = true;
            bool lagrange = true;
            bool mechanics = true;
            Set::Scalar value = 0.0;
            ThresholdType type = ThresholdType::Continuous;
        } threshold;
    } pf;

    struct {
        int on = 0;
        Set::Scalar tstart = 0.0;
        Numeric::Interpolator::Linear<Set::Scalar> mobility;
        Numeric::Interpolator::Linear<Set::Scalar> threshold;
        Set::Field<Set::Scalar> L_mf;
        Set::Field<Set::Scalar> threshold_mf;
    } anisotropic_kinetics;

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
        Set::Scalar tstart = 0.0;
        std::vector<model_type> model;
    } mechanics;
    

};
}
#endif
