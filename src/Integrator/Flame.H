#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <limits>

#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "IC/IC.H"
#include "Set/Set.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "IO/ParmParse.H"
#include "Model/Chamber/Ballistic.H"

#include "Base/Mechanics.H"

#include "Model/Propellant/Constant.H"
#include "Model/Propellant/PowerLaw.H"
#include "Model/Propellant/FullFeedback.H"
#include "Model/Propellant/Homogenize.H"

namespace Integrator
{

class Flame : virtual public Integrator::Integrator,
                virtual public Base::Mechanics<Model::Solid::Finite::NeoHookeanPredeformed>
{

    using model_type = Model::Solid::Finite::NeoHookeanPredeformed;

public:
    static constexpr const char* name = "flame";

    Flame();
    Flame(IO::ParmParse& pp);

    static void Forbids(IO::ParmParse& pp);

    static void Parse(Flame& value, IO::ParmParse& pp);

    virtual ~Flame()
    {
        delete ic_eta;
        delete ic_phi;
        delete thermal.ic_temp;
        delete ic_laser;
        delete bc_eta;
        delete bc_temp;
    }

    // NOTE: members below are public (rather than protected) because their
    // override methods define extended __device__ lambdas, and the parameter
    // structs are captured by value into device lambdas. nvcc forbids both an
    // extended device lambda inside a private/protected member function and a
    // captured variable whose type is a private/protected (or unnamed) member.
    // See docs/gpu_device_capture_conventions.md for the capture convention.
public:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar a_time, int a_iter) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
        const amrex::MFIter& mfi, const amrex::Box& box) override;
    void UpdateModel(int a_step, Set::Scalar a_time) override;

    Set::Field<Set::Scalar> temp_mf;
    Set::Field<Set::Scalar> temp_old_mf;
    Set::Field<Set::Scalar> temps_mf;

    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
    Set::Field<Set::Scalar> psi_mf; // floored copy of eta used to weight the elastic solver
    Set::Field<Set::Scalar> eta_0_mf;
    Set::Field<Set::Scalar> mdot_mf;

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    Set::Field<Set::Scalar> alpha_mf;
    Set::Field<Set::Scalar> heatflux_mf;

    Set::Field<Set::Scalar> laser_mf;

    Set::Field<Set::Scalar> L_mf;

    BC::BC<Set::Scalar>* bc_temp = nullptr;
    BC::BC<Set::Scalar>* bc_eta = nullptr;
    IC::IC<Set::Scalar>* ic_phi = nullptr;
    IC::IC<Set::Scalar>* ic_laser = nullptr;

    Set::Scalar phi_refinement_criterion = std::numeric_limits<Set::Scalar>::infinity();
    Set::Scalar m_refinement_criterion = NAN;
    Set::Scalar t_refinement_criterion = NAN;
    Set::Scalar t_refinement_restriction = NAN;
    Set::Scalar small = NAN;
    IC::IC<Set::Scalar>* ic_eta = nullptr;
    bool plot_field = true;
    int variable_pressure = -1;

    bool homogenized = false;

    amrex::BoxArray prev_finest_ba;
    int prev_finest_level = -1;
    int end_initial_refine = 0;

    // Phase field model variables
    struct pf_struct {
        Set::Scalar eps = NAN;
        Set::Scalar lambda = NAN;
        Set::Scalar kappa = NAN;
        Set::Scalar w1 = NAN, w12 = NAN, w0 = NAN;
        Set::Scalar min_eta = 0.001;
        int relax_steps = 0;
    } pf;

    // Thermal transport variables
    struct thermal_struct {
        bool on = 0;
        Set::Scalar hc = NAN;
        Set::Scalar Tref = NAN;
        Set::Scalar Tfluid = NAN;
        Set::Scalar Tcutoff = NAN;
        Set::Scalar end_initial_refine_time = NAN;
        Set::Scalar phi_refinement_criterion_inital = NAN;
        Set::Field<Set::Scalar> has_exceeded_Tcutoff;
        IC::IC<Set::Scalar>* ic_temp = nullptr;
    } thermal;

    // Thermoelastic variables
    struct elastic_struct {
        int on = 0;
        Set::Scalar Telastic = NAN;
        // Grain-scale (homogenized propellant) elastic models: a single
        // homogenized propellant material plus the void (combustion gas) and
        // casing (confinement) phases for the three-material partition.
        model_type model_prop, model_void, model_casing;
        // Mesoscale (resolved AP/HTPB propellant) elastic models, blended by the
        // species field phi. Used instead of model_prop when the propellant is
        // not homogenized (e.g. propellant.type = fullfeedback). See UpdateModel.
        model_type model_ap, model_htpb;
        Set::Scalar traction = NAN;
        // When 1, the interface traction follows the evolving chamber pressure
        // instead of the constant elastic.traction value.
        int traction_from_chamber = 0;
        int phirefinement = -1;
        // Floor for the psi weighting field used by the elastic solver. psi=eta
        // zeroes the gas region, which leaves the masked operator near-singular on
        // thin/complex solid features and makes MLMG diverge. With a floor f, the
        // solver sees psi = f + (1-f)*eta, keeping the operator positive-definite in
        // the gas (whose stress we ignore). Default 0 reproduces the old psi=eta.
        // Only used when use_psi=1.
        Set::Scalar psi_floor = 0.0;
        // Whether to mask the elastic operator with the psi field at all. When 0,
        // no psi is installed (psi_avg=1 everywhere) and the gas is conditioned
        // purely by its real, unitized model_void stiffness -- a cleaner model with
        // no arbitrary scalar multiplying the void modulus. The masked operator with
        // (psi_floor*kappa_void) is numerically equivalent to use_psi=0 with
        // kappa_void scaled to that product, so a soft physical void replaces both.
        int use_psi = 1;
    } elastic;

    // Chamber state variables (not inputs)
    struct chamber_struct {
        // Integrated variables
        Set::Scalar volume = 0.0;
        Set::Scalar area = 0.0;
        //Set::Scalar massflux = 0.0;
        Set::Scalar mdot = 0.0;
        // Intensive variables
        Set::Scalar pressure = NAN;

        Model::Chamber::Ballistic  model;

    } chamber;

    // Thermo output diagnostics (filled in TimeStepComplete, written to thermo.dat)
    Set::Scalar thermo_max_temp     = 0.0;
    Set::Scalar thermo_mdot_max     = 0.0;
    Set::Scalar thermo_heatflux_max = 0.0;
    Set::Scalar thermo_L_max        = 0.0;
    Set::Scalar thermo_eta_min      = 1.0;

    // Propellant model
    Model::Propellant::Propellant<
        Model::Propellant::Constant,
        Model::Propellant::PowerLaw,
        Model::Propellant::FullFeedback,
        Model::Propellant::Homogenize> propellant;

};
}

#endif
