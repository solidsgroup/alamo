#ifndef INTEGRATOR_NARROWBANDLEVELSET_H // Include guards
#define INTEGRATOR_NARROWBANDLEVELSET_H // 

// Standard library includes
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>
#include <stdexcept>
#include <memory>

// Alamo includes
#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IC/IC.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"

// Numeric includes
//#include "Numeric/FluxHandler.H"
//#include "Numeric/TimeStepper.H"
#include "Numeric/Stencil.H"

namespace Numeric {
    template <typename T> class FluxReconstruction;
    template <typename T> class FluxMethod;
    template <typename T> class FluxHandler;  
    template <typename T> class FirstOrderReconstruction;
    template <typename T> class LocalLaxFriedrichsMethod;
    template <typename T> class TimeStepper;
    template <typename T> class TimeSteppingScheme;
    template <typename T> class EulerForwardScheme;  
    template <typename T> class RK3Scheme;
}

namespace Integrator
{
class NarrowBandLevelset : virtual public Integrator
{
public:

    // Constructors
    NarrowBandLevelset(int a_nghost = 2); // Default constructor
    NarrowBandLevelset(IO::ParmParse& pp); // Constructor that triggers parse
    static void Parse(NarrowBandLevelset& value, IO::ParmParse& pp);

    // Destructor
    virtual ~NarrowBandLevelset()
    {
        delete ic_ls;
        delete ic_velocity;
        
        delete bc_ls;
        
    }
    
    // Public getters to access protected/private members in Numeric templates
    int GetNumGhostCells() const { return number_of_ghost_cells; }
    const Set::Field<Set::Scalar>& GetLevelSetField() const { return ls_mf; }
    const Set::Field<Set::Scalar>& GetLevelSetFieldOld() const { return ls_old_mf; }
    const Set::Field<Set::Scalar>& GetVelocityField() const { return velocity_mf; }
    const Set::Field<Set::Scalar>& GetXFluxField() const { return XFlux_mf; }
#if AMREX_SPACEDIM >= 2
    const Set::Field<Set::Scalar>& GetYFluxField() const { return YFlux_mf; }
#endif
#if AMREX_SPACEDIM == 3
    const Set::Field<Set::Scalar>& GetZFluxField() const { return ZFlux_mf; }
#endif
    
    // Flux Reconstruction Schemes
    enum class FluxReconstruction {
        FirstOrder,               // Basic first-order reconstruction
        ThirdOrderENO,            // Third-order Essentially Non-Oscillatory (ENO)
        FifthOrderWENO,           // Fifth-order Weighted Essentially Non-Oscillatory (WENO)
        FifthOrderWENOWithChar    // Fifth-order WENO with Characteristic decomposition
    };

    // Flux Schemes
    enum class FluxScheme {
        LocalLaxFriedrichs,   // Local Lax-Friedrichs (Rusanov) flux
        HLLCFlux             // Harten-Lax-van Leer Contact (HLLC) flux
    };

    enum class TimeSteppingScheme {
        ForwardEuler,           // First-order explicit Euler method
        RK3                    // Third-order Runge-Kutta method
    };

    // Subclass to manage feature maps
    class FeatureMaps {
    public:
        // Singleton access to ensure one instance
        static FeatureMaps& getInstance() {
            static FeatureMaps instance;
            return instance;
        }

        // Accessor methods for the maps
        const std::map<std::string, FluxReconstruction>& getFluxReconstructionMap() const {
            return fluxReconstructionMap;
        }

        const std::map<std::string, FluxScheme>& getFluxSchemeMap() const {
            return fluxSchemeMap;
        }

        const std::map<std::string, TimeSteppingScheme>& getTimeSteppingSchemeMap() const {
            return timeSteppingSchemeMap;
        }

    private:
        // Private constructor to initialize the maps
        FeatureMaps() {
            fluxReconstructionMap = {
                {"FirstOrder", FluxReconstruction::FirstOrder},
                {"ThirdOrderENO", FluxReconstruction::ThirdOrderENO},
                {"FifthOrderWENO", FluxReconstruction::FifthOrderWENO},
                {"FifthOrderWENOWithChar", FluxReconstruction::FifthOrderWENOWithChar}
            };

            fluxSchemeMap = {
                {"LocalLaxFriedrichs", FluxScheme::LocalLaxFriedrichs},
                {"HLLCFlux", FluxScheme::HLLCFlux}
            };

            timeSteppingSchemeMap = {
                {"ForwardEuler", TimeSteppingScheme::ForwardEuler},
                {"RK3", TimeSteppingScheme::RK3}
            };
        }

        // Disable copy and assignment
        FeatureMaps(const FeatureMaps&) = delete;
        FeatureMaps& operator=(const FeatureMaps&) = delete;

        // The maps for different schemes
        std::map<std::string, FluxReconstruction> fluxReconstructionMap;
        std::map<std::string, FluxScheme> fluxSchemeMap;
        std::map<std::string, TimeSteppingScheme> timeSteppingSchemeMap;
    };

    // Accessors to the singleton
    static const FeatureMaps& getFeatureMaps() {
        return FeatureMaps::getInstance();
    }
    
protected:
    // Variables (multi-fabs)
    Set::Field<Set::Scalar>    ls_mf;                   // Levelset field variable (current timestep)
    Set::Field<Set::Scalar>    ls_old_mf;               // Levelset field variable (previous timestep)
    Set::Field<Set::Scalar>    velocity_mf;             // Velocity field variable (Assumed constant for now)
    Set::Field<Set::Scalar>    XFlux_mf;                // X-direction flux field
    Set::Field<Set::Scalar>    YFlux_mf;                // Y-direction flux field
    Set::Field<Set::Scalar>    ZFlux_mf;                // Z-direction flux field

    // Solver setup and control functions (override)
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;
    
    // Additional functions
    void Advect(int lev, Set::Scalar time, Set::Scalar dt);
    void UpdateInterfaceVelocity(int lev);
    void Reinitialize(int lev, Set::Scalar time);
    void ApplyNarrowBanding(int lev);
    void ComputeAndSetNewTimeStep();
    Set::Scalar GetTimeStep();
    void ApplyBoundaryConditions(int lev, Set::Scalar time);

    std::shared_ptr<Numeric::FluxHandler<NarrowBandLevelset>> fluxHandler;
    std::shared_ptr<Numeric::TimeStepper<NarrowBandLevelset>> timeStepper;

    template <typename T>
    friend class Numeric::FluxReconstruction;

    template <typename T>
    friend class Numeric::FluxMethod;    

    template <typename T>
    friend class Numeric::FluxHandler;    

    template <typename T>
    friend class Numeric::FirstOrderReconstruction;

    template <typename T>
    friend class Numeric::LocalLaxFriedrichsMethod;

    template <typename T>
    friend class Numeric::TimeSteppingScheme;

    template <typename T>
    friend class Numeric::EulerForwardScheme;

    template <typename T>
    friend class Numeric::RK3Scheme;    

private:
    //
    // Definition of parameters set only at instantiation by
    // constructors. 
    //
    int current_timestep            = 0;
    const int number_of_components  = 1;     // Number of components (1-levelset value)
    const int number_of_ghost_cells = 2;     // Number of ghost cells

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC<Set::Scalar>* ic_ls = nullptr;                     // Object used to initialize levelset field
    IC::IC<Set::Scalar>* ic_velocity = nullptr;               // Object used to initialize velocity field
    BC::BC<Set::Scalar>* bc_ls = nullptr;        // Object used to update temp levelset field boundary ghost cells
    BC::BC<Set::Scalar>* bc_velocity = nullptr;  // Object used to update temp velocity field boundary ghost cells
    BC::Nothing bc_nothing;
    
    FluxReconstruction reconstruction_method;    // Selected flux reconstruction scheme
    FluxScheme flux_scheme;                      // Selected flux scheme
    TimeSteppingScheme temporal_scheme;          // Selected time-stepping scheme
    
    Set::Scalar cflNumber;                       // CFL number (read from input)

}; //NarrowBandLevelset Class

} // namespace Integrator
#endif
