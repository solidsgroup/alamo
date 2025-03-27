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
        delete bc_ls;
        
    }
    
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
    // Each level set has its own material ID and phi field
    struct LevelSetData
    {
        int id = -1;
        int material_id = -1;
        Set::Field<Set::Vector> velocity_mf;
        Set::Field<Set::Vector> normal_mf;
        Set::Field<Set::Scalar> curvature_mf;
        
        // KEY CONCEPT: Each level set tracks its own narrow band region
        // as a collection of boxes - enabling quick intersection tests
        bool has_narrowband;
        amrex::BoxArray narrowband_boxes;
        amrex::DistributionMapping narrowband_dm;
    };
    
    // Variables (multi-fabs)
    Set::Field<Set::Scalar>    ls_mf;                   // Levelset field variable (current timestep)
    Set::Field<Set::Scalar>    ls_old_mf;               // Levelset field variable (previous timestep)
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

    // Class specific functions
    void ComputeNarrowBandBox(int lev, int ls_id);
    void ComputeNarrowBandMapping(int ls_id);
    void UpdateNarrowband(int lev, int ls_id);
    bool BoxIntersectsNarrowBand(const amrex::Box& box, int ls_id, int lev) const;
    void ComputeNormal(int lev, int ls_id);
    void ComputeCurvature(int lev, int ls_id);
    void ComputeGeometryQuantities(int lev, int ls_id);
    void ComputeAndSetNewTimeStep();
    Set::Scalar GetTimeStep();
    void Advect(int lev, Set::Scalar time, Set::Scalar dt);
    void UpdateInterfaceVelocity(int lev);
    void Reinitialize(int lev);

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
    
    // Define LevelSetData structure object
    std::vector<LevelSetData> level_sets;
    
    // Define size for narrowband
    Set::Scalar narrow_band_width = 6.0;

    // Define constant velocity array
    Set::Vector constant_velocity = Set::Vector::Zero();

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC<Set::Scalar>* ic_ls = nullptr;        // Object used to initialize levelset field
    BC::BC<Set::Scalar>* bc_ls = nullptr;        // Object used to update temp levelset field boundary ghost cells
    BC::Nothing bc_nothing;
    
    FluxReconstruction reconstruction_method;    // Selected flux reconstruction scheme
    FluxScheme flux_scheme;                      // Selected flux scheme
    TimeSteppingScheme temporal_scheme;          // Selected time-stepping scheme

    Set::Scalar cflNumber = 0.5;                 // CFL number (read from input)

}; //NarrowBandLevelset Class

} // namespace Integrator
#endif
