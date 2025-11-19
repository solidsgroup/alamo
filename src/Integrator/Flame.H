#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <limits>
#include <vector>

#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "IC/IC.H"
#include "Set/Set.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "IO/ParmParse.H"

#include "Base/Mechanics.H"

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
        delete thermal.ic_temp;
        delete ic_laser;
        delete bc_temp;
    }

protected:
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override;
    void TimeStepComplete(Set::Scalar a_time, int a_iter) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;
    void Regrid(int lev, Set::Scalar time) override;
    void Integrate(int amrlev, Set::Scalar time, int step,
        const amrex::MFIter& mfi, const amrex::Box& box) override;
    void UpdateModel(int a_step, Set::Scalar a_time) override;

private:
    int num_species = 0;

    Set::Field<Set::Scalar> temp_mf;
    Set::Field<Set::Scalar> temp_old_mf;
    Set::Field<Set::Scalar> temps_mf;

    Set::Field<Set::Scalar> mdot_mf;

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> phi_old_mf;
    Set::Field<Set::Scalar> field;
    Set::Field<Set::Scalar> alpha_mf;
    Set::Field<Set::Scalar> heatflux_mf;

    Set::Field<Set::Scalar> laser_mf;

    std::vector<BC::BC<Set::Scalar>*> bc_phi;
    std::vector<IC::IC<Set::Scalar>*> ic_phi;
    BC::BC<Set::Scalar>* bc_temp = nullptr;
    IC::IC<Set::Scalar> *ic_laser = nullptr;

    Set::Scalar phi_refinement_criterion = std::numeric_limits<Set::Scalar>::infinity();
    Set::Scalar m_refinement_criterion = NAN;
    Set::Scalar t_refinement_criterion = NAN;
    Set::Scalar t_refinement_restriction = NAN;
    Set::Scalar small = NAN;
    bool plot_field = true;
    int variable_pressure = -1;

    // Phase field model variables
    struct {
        Set::Scalar eps = NAN;
        Set::Scalar lambda = NAN;
        Set::Scalar kappa = NAN;
        Set::Scalar w1 = NAN, w12 = NAN, w0 = NAN;
        Set::Scalar min_eta = 0.001;  // update name?
    } pf;

    // Thermal transport variables
    struct {
        bool on = 0;
        Set::Scalar hc = NAN;
        Set::Scalar Tref = NAN;
        Set::Scalar Tfluid = NAN;
        IC::IC<Set::Scalar>* ic_temp = nullptr;
    } thermal;

    // Thermoelastic variables
    struct {
        int on = 0;
        Set::Scalar Telastic = NAN;
        model_type model_ap, model_htpb;
        Set::Scalar traction = NAN;
        int phirefinement = -1;  // update name?
    } elastic;

    // Chamber state variables (not inputs)
    struct {
        // Integrated variables
        Set::Scalar volume = 0.0;
        Set::Scalar area = 0.0;
        Set::Scalar massflux = 0.0;
        // Intensive variables
        Set::Scalar pressure = NAN;
    } chamber;


    // Propellant model
    Model::Propellant::Propellant<
        Model::Propellant::PowerLaw,
        Model::Propellant::FullFeedback,
        Model::Propellant::Homogenize> propellant;
};
}

#endif
