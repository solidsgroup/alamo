#ifndef INTEGRATOR_LOWMACH_H
#define INTEGRATOR_LOWMACH_H

#include "Integrator/Integrator.H"
#include "IO/ParmParse.H"

#include "BC/BC.H"
#include "BC/Constant.H"
#include "BC/Expression.H"
#include "BC/Nothing.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/IC.H"
#include "IC/PNG.H"
#include "IC/PSRead.H"
#include "Model/Gas/Gas.H"
#include "Model/Chemistry/Chemistry.H"
#include "Model/Chemistry/FiniteRate.H"
#include "Model/Chemistry/Frozen.H"
#include "Model/Chemistry/GrossModel.H"
#include "Model/Chemistry/GrossModel_Aluminized.H"
#include "Model/Mechanism/Mechanism.H"
#include "Model/Mechanism/PhaseChange.H"
#include "Model/Mechanism/InterphaseReaction.H"
#include "Model/Mechanism/InterfacialHeatSource.H"
#include "Model/PhaseField/MultiphaseInterface.H"
#include "Model/Solid/Finite/NeoHookean.H"
#include "Set/Base.H"
#include "Numeric/Advect/Advect.H"
#include "Numeric/Advect/Centered.H"
#include "Numeric/Advect/MUSCL.H"
#include "Numeric/Advect/QUICK.H"
#include "Numeric/Advect/Upwind.H"
#include "Numeric/Advect/WENO5.H"
#include "Numeric/ReferenceMap/Reconstruction.H"
#include "Operator/Diffusion.H"
#include "Operator/PressurePoisson.H"

namespace Integrator
{
class LowMach : virtual public Integrator
{
public:
    static constexpr const char* name = "lowmach";
    using Chemistry = Model::Chemistry::Chemistry<
        Model::Chemistry::Frozen,
        Model::Chemistry::FiniteRate,
        Model::Chemistry::GrossModel,
        Model::Chemistry::GrossModel_Aluminized>;

    LowMach() {};
    LowMach(IO::ParmParse& pp);
    static void Parse(LowMach& value, IO::ParmParse& pp);

    virtual ~LowMach()
    {
        delete velocity_bc;
        delete temperature_bc;
        delete component_density_bc;
        delete pressure_bc;
        delete xi_bc;

        delete velocity_ic;
        delete temperature_ic;
        delete heat_source_ic;
        for (auto* ic : component_density_ic)
            delete ic;
        delete pressure_ic;
        delete xi_ic;
    }

public:
    // NVCC requires enclosing functions for extended device lambdas to be public.
    void Initialize(int lev) override;
    void Regrid(int lev, Set::Scalar time) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int iter) override;
    void PreparePlotFile(Set::Scalar time, const amrex::Vector<int>& iter) override;
    void PrintDiagnostics(Set::Scalar time, int iter);
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

    void UpdateSolidStress(int lev,
                            const amrex::MultiFab& u_mf,
                            const amrex::MultiFab& eta_mf,
                            const amrex::MultiFab& xi_mf,
                            bool write_diagnostics = false);
    void UpdateComponentState(int lev, const amrex::MultiFab& component_density_mf);
    void UpdateInterfacialChemicalPotential(int lev);
    void ApplyInterfacialTransport(Set::Scalar time, Set::Scalar dt);
    void UpdateDerivedDiagnostics(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf);
    void AdvanceChemistry(int lev, amrex::MultiFab& temperature_mf,
                            amrex::MultiFab& component_density_mf, Set::Scalar dt);
    void ApplyImplicitPhaseChange(Set::Scalar time, Set::Scalar dt);
    void ApplyEquilibriumPhaseChange(Set::Scalar time, Set::Scalar dt);
    void ApplyKineticPhaseChange(Set::Scalar time, Set::Scalar dt);
    void ApplyImplicitDiffusion(Set::Scalar time, Set::Scalar dt);
    void UpdateFreeRigidBodyStates();

    using ThermalData = amrex::GpuTuple<
        Model::Gas::Gas::DeviceData, int, int, Set::Scalar, Set::Scalar, bool,
        const Set::Scalar*, const Set::Scalar*, const Set::Scalar*,
        const Set::Scalar*, const Set::Scalar*>;
    using ChemistryData = amrex::GpuTuple<
        bool, Model::Chemistry::GrossModel, Chemistry::ODESolver,
        const Chemistry*, const Model::Gas::Gas*>;
    using ThermochemicalData = amrex::GpuTuple<
        ThermalData, ChemistryData, bool, bool, bool>;

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static amrex::GpuTuple<
        Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar>
    ComputeThermalState(
        Set::Patch<const Set::Scalar> component_density,
        Set::Scalar temperature, int i, int j, int k,
        const ThermalData& data);

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static Set::Scalar ComputeViscosity(
        Set::Patch<const Set::Scalar> component_density,
        Set::Scalar temperature, int i, int j, int k,
        const ThermalData& data);

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static
    amrex::GpuTuple<Model::Chemistry::SpeciesArray, Set::Scalar /*temperature*/, Set::Scalar /*dilatation*/>
    ComputeThermochemicalSource(
        Set::Patch<const Set::Scalar> component_density,
        Set::Patch<const Set::Scalar> temperature,
        int i, int j, int k, const Set::Scalar* dx, Set::Scalar dt,
        const ThermochemicalData& data,
        bool include_reaction = true,
        bool project_negative_gas_density = false);
    void ProjectVelocity(Set::Scalar time, Set::Scalar dt);
    void RHS(int lev, Set::Scalar time, Set::Scalar dt,
            amrex::MultiFab& u_rhs_mf,
            amrex::MultiFab& T_rhs_mf,
            amrex::MultiFab& component_density_rhs_mf,
            amrex::MultiFab* xi_rhs_mf,
            const amrex::MultiFab& u_mf,
            const amrex::MultiFab& T_mf,
            const amrex::MultiFab& component_density_mf,
            const amrex::MultiFab* xi_mf);

protected:
    Set::Field<Set::Scalar> velocity_mf;
    Set::Field<Set::Scalar> velocity_old_mf;
    Set::Field<Set::Scalar> temperature_mf;
    Set::Field<Set::Scalar> temperature_old_mf;
    Set::Field<Set::Scalar> heat_source_mf;
    Set::Field<Set::Scalar> component_density_mf;
    Set::Field<Set::Scalar> component_density_old_mf;
    Set::Field<Set::Scalar> mass_fraction_mf;
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> rigid_eta_mf;
    Set::Field<Set::Scalar> rigid_species_eta_mf;
    Set::Field<Set::Scalar> fixed_rigid_eta_mf;
    Set::Field<Set::Scalar> liquid_species_eta_mf;
    Set::Field<Set::Scalar> interfacial_volume_fraction_mf;
    Set::Field<Set::Scalar> gas_volume_fraction_mf;
    Set::Field<Set::Scalar> interfacial_chemical_potential_mf;
    Set::Field<Set::Scalar> xi_mf;
    Set::Field<Set::Scalar> xi_old_mf;

    Set::Field<Set::Scalar> density_mf;
    Set::Field<Set::Scalar> momentum_mf;
    Set::Field<Set::Scalar> pressure_mf;
    Set::Field<Set::Scalar> pressure_correction_mf;
    Set::Field<Set::Scalar> chemistry_dilatation_mf;
    Set::Field<Set::Scalar> phase_change_dilatation_mf;
    Set::Field<Set::Scalar> phase_change_heat_mf;
    Set::Field<Set::Scalar> kinetic_diffuse_moment_mf;
    Set::Field<Set::Scalar> diffusion_dilatation_mf;
    Set::Field<Set::Scalar> energy_mf;
    Set::Field<Set::Scalar> mole_fraction_mf;
    Set::Field<Set::Scalar> vorticity_mf;
    Set::Field<Set::Scalar> viscosity_mf;
    Set::Field<Set::Scalar> thermal_conductivity_coeff_mf;
    Set::Field<Set::Scalar> diffusion_coeff_mf;
    Set::Field<Set::Scalar> wdot_mf;
    Set::Field<Set::Scalar> qdot_mf;
    Set::Field<Set::Matrix> deformation_gradient_mf;
    Set::Field<Set::Matrix> solid_deviatoric_stress_mf;

    BC::BC<Set::Scalar>* velocity_bc = nullptr;
    BC::BC<Set::Scalar>* temperature_bc = nullptr;
    BC::BC<Set::Scalar>* component_density_bc = nullptr;
    BC::BC<Set::Scalar>* pressure_bc = nullptr;
    BC::BC<Set::Scalar>* xi_bc = nullptr;
    BC::Nothing bc_nothing;

    IC::IC<Set::Scalar>* velocity_ic = nullptr;
    IC::IC<Set::Scalar>* temperature_ic = nullptr;
    IC::IC<Set::Scalar>* heat_source_ic = nullptr;
    std::vector<IC::IC<Set::Scalar>*> component_density_ic;
    IC::IC<Set::Scalar>* pressure_ic = nullptr;
    IC::IC<Set::Scalar>* xi_ic = nullptr;

    Model::Gas::Gas gas;
    amrex::Gpu::DeviceVector<Set::Scalar> gas_device_storage;
    amrex::Gpu::DeviceVector<Set::Scalar> thermal_device_storage;
    Model::Gas::Gas::DeviceData gas_device_data;
    ThermalData thermal_data;
    Chemistry chemistry;
    ChemistryData chemistry_device_data;
    ThermochemicalData thermochemical_data;
    Model::Solid::Finite::NeoHookean finite_solid_model;
    using Mechanism = Model::Mechanism::Mechanism<
        Model::Mechanism::PhaseChange,
        Model::Mechanism::InterphaseReaction,
        Model::Mechanism::InterfacialHeatSource>;
    std::vector<Mechanism> mechanisms;
    std::vector<std::string> species_names;
    std::vector<Set::Scalar> reference_density;
    Model::Chemistry::SpeciesArray condensed_specific_heat;
    Model::Chemistry::SpeciesArray condensed_thermal_conductivity;
    Model::Chemistry::SpeciesArray condensed_inverse_reference_density;
    Model::Chemistry::SpeciesArray condensed_temperature_override;
    Model::Chemistry::SpeciesArray condensed_dynamic_viscosity;
    Model::Chemistry::SpeciesArray liquid_inverse_reference_density;
    int ngas_species = 0;
    int nspecies = 0;
    int deformable_solid_species = -1;
    std::vector<int> rigid_solid_species;
    std::vector<Set::Scalar> rigid_relaxation_time;
    std::vector<Set::Vector> rigid_velocity;
    std::vector<bool> rigid_anderson_coupling;
    std::vector<int> rigid_coupling_max_iterations;
    std::vector<Set::Scalar> rigid_coupling_relative_tolerance;
    std::vector<Set::Scalar> rigid_coupling_absolute_tolerance;
    std::vector<int> fixed_rigid_solid_species;
    std::vector<int> free_rigid_solid_species;
    std::vector<int> free_rigid_solid_components;
    std::vector<int> liquid_species;
    std::vector<int> interfacial_solid_species;
    Model::PhaseField::MultiphaseInterface interfacial_model;
    enum class InterfacialTransport
    {
        Advected,
        ConservativeAllenCahn
    };
    InterfacialTransport interfacial_transport =
        InterfacialTransport::Advected;

    Set::Vector g = Set::Vector::Zero();
    Set::Scalar cfl = NAN;
    Set::Scalar phase_field_cfl = NAN;
    enum class ChemistryTimestepMode
    {
        Off,
        Report,
        Limit
    };
    ChemistryTimestepMode chemistry_timestep_mode =
        ChemistryTimestepMode::Off;
    // Both controls below are normalized fractions and therefore
    // intentionally nondimensional.
    Set::Scalar chemistry_max_fractional_change = 0.1;
    Set::Scalar chemistry_reactant_mass_fraction_floor = 1.0e-4;
    Set::Scalar chemistry_timescale_min = NAN;
    Set::Scalar chemistry_timestep_candidate = NAN;
    Set::Scalar cfl_v = 1.0e100;
    Set::Scalar small = 1.0e-12;
    Set::Scalar density_floor = 1.0e-12;
    Set::Scalar pressure_floor = 1.0e-12;
    Set::Scalar pressure_scale = 1.0;
    Set::Scalar pressure_reference = NAN;
    Set::Scalar velocity_refinement_criterion = NAN;
    Set::Scalar pressure_refinement_criterion = NAN;
    Set::Scalar temperature_refinement_criterion = NAN;
    Set::Scalar reaction_refinement_criterion = NAN;
    Set::Scalar eta_refinement_criterion = NAN;
    Set::Scalar reinitialize_condensed_composition_eta_min = 0.99;
    int diagnostics_interval = 0;
    bool diagnostics_extended_fields = false;
    bool reinitialize_condensed_composition = false;
    bool projection_enabled = true;
    bool projection_update_pressure = true;
    bool include_viscosity = true;
    bool include_conduction = true;
    // Enabled by default when liquid species are configured; it must remain
    // false for gas/solid-only systems because no capillary fields exist.
    bool interfacial_forces_enabled = false;
    bool liquid_solid_interface_enabled = false;
    bool implicit_momentum_diffusion = false;
    bool implicit_thermal_diffusion = false;
    bool implicit_species_diffusion = false;
    bool common_species_diffusivity = false;
    bool has_equilibrium_phase_change = false;
    bool has_kinetic_phase_change = false;
    bool has_split_phase_change = false;
    bool condensed_thermal_transport = false;
    bool advect_temperature = true;
    Set::Scalar finite_solid_reference_density = 1.0;
    Set::Scalar finite_solid_eta_threshold = 0.5;
    Set::Scalar finite_solid_J_floor = 1.0e-6;
    Set::Scalar finite_solid_viscosity = 0.0;
    Set::Scalar finite_solid_interface_viscosity = 0.0;
    Set::Scalar finite_solid_deviatoric_stress_divergence_sign = 0.0;
    Set::Scalar interfacial_reference_density_min = 1.0e100;
    Set::Scalar interfacial_thickness = 0.0;
    Set::Scalar conservative_allen_cahn_diffusivity = 0.0;
    Set::Scalar previous_free_rigid_dt = NAN;

    struct FreeRigidBodyState
    {
        Set::Scalar mass = 0.0;
        Set::Vector center = Set::Vector::Zero();
        Set::Vector velocity = Set::Vector::Zero();
        Set::Vector angular_velocity = Set::Vector::Zero();
        Set::Matrix inertia = Set::Matrix::Zero();
        Set::Scalar radius_of_gyration = 0.0;
        bool valid = false;
    };
    std::vector<int> last_free_rigid_coupling_iterations;
    std::vector<Set::Scalar> last_free_rigid_coupling_residual;
    std::vector<bool> last_free_rigid_coupling_converged;
    std::vector<Set::Scalar> free_rigid_body_time;
    std::vector<FreeRigidBodyState> free_rigid_bodies;
    std::vector<FreeRigidBodyState> previous_free_rigid_projection_increments;

    Numeric::Advect::Advect<Numeric::Advect::MUSCL,
                            Numeric::Advect::Upwind,
                            Numeric::Advect::Centered,
                            Numeric::Advect::QUICK,
                            Numeric::Advect::WENO5> advect;
    Numeric::ReferenceMap::Reconstruction reference_map_reconstruction;
    Operator::Diffusion diffusion;
    Operator::PressurePoisson pressure_poisson;
    bool synchronize_restart_state = false;
};
}

#endif
