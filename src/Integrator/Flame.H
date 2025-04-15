#ifndef INTEGRATOR_FLAME_H
#define INTEGRATOR_FLAME_H

#include <limits>

#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "IC/IC.H"
#include "Set/Set.H"
#include "Model/Solid/Finite/NeoHookeanPredeformed.H"
#include "IO/ParmParse.H"

#include "Base/Mechanics.H"

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

    Set::Field<Set::Scalar> temp_mf;
    Set::Field<Set::Scalar> temp_old_mf;
    Set::Field<Set::Scalar> temps_mf;
    Set::Field<Set::Scalar> temps_old_mf;

    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
    Set::Field<Set::Scalar> mob_mf;
    Set::Field<Set::Scalar> mdot_mf;

    Set::Field<Set::Scalar> phi_mf;
    Set::Field<Set::Scalar> field;
    Set::Field<Set::Scalar> alpha_mf;
    Set::Field<Set::Scalar> heatflux_mf;

    Set::Field<Set::Scalar> laser_mf;

    BC::BC<Set::Scalar>* bc_temp = nullptr;
    BC::BC<Set::Scalar>* bc_eta = nullptr;
    IC::IC<Set::Scalar>* ic_phi = nullptr;
    IC::IC<Set::Scalar>* ic_laser = nullptr;

    Set::Scalar phi_refinement_criterion = std::numeric_limits<Set::Scalar>::infinity();
    Set::Scalar m_refinement_criterion = NAN;
    Set::Scalar t_refinement_criterion = NAN;
    Set::Scalar t_refinement_restriction = NAN;
    Set::Scalar zeta = NAN;
    Set::Scalar zeta_0 =NAN;
    Set::Scalar small = NAN;
    Set::Scalar base_time = NAN;
    IC::IC<Set::Scalar>* ic_eta = nullptr;
    int ghost_count = -1;
    int homogeneousSystem = -1;
    bool plot_field = true;
    int variable_pressure = -1;

    struct {
        Set::Scalar gamma = NAN;
        Set::Scalar eps = NAN;
        Set::Scalar lambda = NAN;
        Set::Scalar kappa = NAN;
        Set::Scalar w1 = NAN, w12 = NAN, w0 = NAN;
        Set::Scalar min_eta = 0.001;
        int evolve = 1;
    } pf;

    struct {
        Set::Scalar P = NAN;
        Set::Scalar gamma = NAN;

        struct {
            Set::Scalar a1 = NAN, a2 = NAN, a3 = NAN;
            Set::Scalar b1 = NAN, b2 = NAN, b3 = NAN;
            Set::Scalar c1 = NAN;
            int mob_ap = NAN;
            int dependency = 1;
            Set::Scalar h1 = NAN, h2 = NAN;

        } arrhenius;
        struct {
            Set::Scalar r_ap = NAN, r_htpb = NAN, r_comb = NAN;
            Set::Scalar n_ap = NAN, n_htpb = NAN, n_comb = NAN;
            Set::Scalar a_fit = NAN, b_fit = NAN, c_fit = NAN;
        } power;
    } pressure;

    struct {
        bool on = 0;
        Set::Scalar rho_ap = NAN, rho_htpb = NAN;
        Set::Scalar k_ap = NAN, k_htpb = NAN;
        Set::Scalar cp_ap = NAN, cp_htpb = NAN;
        Set::Scalar q0 = NAN;
        Set::Scalar bound = NAN;
        Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
        Set::Scalar E_ap = NAN, E_htpb = NAN;
        Set::Scalar hc = NAN, w1 = NAN;
        Set::Scalar modeling_ap = NAN, modeling_htpb = NAN;
        Set::Scalar mlocal_ap = NAN, mlocal_htpb = NAN, mlocal_comb = NAN;
        Set::Scalar T_fluid = NAN;
        Set::Scalar massfraction = NAN;
        Set::Scalar disperssion1 = NAN;
        Set::Scalar disperssion2 = NAN;
        Set::Scalar disperssion3 = NAN;
        IC::IC<Set::Scalar>* ic_temp = nullptr;
    } thermal;

    struct {
        int on = 0;
        Set::Scalar Tref = NAN;
        model_type model_ap, model_htpb;
        Set::Scalar traction = NAN;
        int phirefinement = -1;
    } elastic;

    // Integrated variables
    Set::Scalar volume = 0.0;
    Set::Scalar area = 0.0;
    Set::Scalar chamber_area = 0.0;
    Set::Scalar chamber_pressure = 1.0;
    Set::Scalar massflux = 0.00000001;

    Set::Scalar x_len = NAN;
    Set::Scalar y_len = NAN;

    BC::BC<Set::Scalar>* bc_psi = nullptr;
};
}

#endif
