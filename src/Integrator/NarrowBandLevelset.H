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
    NarrowBandLevelset(int a_nghost = 1); // Default constructor
    NarrowBandLevelset(IO::ParmParse& pp); // Constructor that triggers parse
    static void Parse(NarrowBandLevelset& value, IO::ParmParse& pp);

    // Destructor
    virtual ~NarrowBandLevelset()
    {
        delete ic_ls;
        delete ic_velocity;
        delete bc_ls;
        delete bc_velocity;
        
    }

    // Getter for Narrowband tube type
    AMREX_FORCE_INLINE static int GetOutsideNarrowBandPos() {
        return NarrowBandTubeType::OutsideNarrowBandPos;
    }

    AMREX_FORCE_INLINE static int GetOuterEdge() {
        return NarrowBandTubeType::OuterEdge;
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
        Set::Field<Set::Scalar> velocity_mf;
        Set::Field<Set::Scalar> normal_mf;
        Set::Field<Set::Scalar> curvature_mf;
        
        // KEY CONCEPT: Each level set tracks its own narrow band region
        // as a collection of boxes - enabling quick intersection tests
        bool has_narrowband;

        std::unique_ptr<amrex::iMultiFab> Tube_imf;
        //amrex::Vector<amrex::IntVect> narrowband_cells; 
        //amrex::BoxList narrowband_boxes;
        amrex::BoxArray narrowband_ba;
        amrex::DistributionMapping narrowband_dm;
        amrex::Vector<int> narrowband_flags;
    };

    // Define a structure for the values of the narrowband
    struct NarrowBandTubeType
    {
        static constexpr int Interface =  0;            // 0 level set
        static constexpr int InnerTube = -1;            // phi < 0 && phi >= -6.0 * min_DX
        static constexpr int OuterTube =  1;            // phi > 0 && phi <= 6.0 * min_DX
        static constexpr int InnerEdge = -2;            // neighbors of narrowband (inside)
        static constexpr int OuterEdge =  2;            // neighbors of narrowband (outside)
        static constexpr int OutsideNarrowBandNeg = -3; // phi > 6.0 * min_DX
        static constexpr int OutsideNarrowBandPos = 3;  // Phi < -6.0 * min_DX
    };

    struct Neighbors {
        static constexpr int num_neighbors =
        #if AMREX_SPACEDIM == 1
            2;
        #elif AMREX_SPACEDIM == 2
            4;
        #elif AMREX_SPACEDIM == 3
            6;
        #endif
    
        static constexpr int offsets[num_neighbors][3] =
        #if AMREX_SPACEDIM == 1
        {
            {-1, 0, 0},
            { 1, 0, 0}
        };
        #elif AMREX_SPACEDIM == 2
        {
            {-1, 0, 0},
            { 1, 0, 0},
            { 0,-1, 0},
            { 0, 1, 0}
        };
        #elif AMREX_SPACEDIM == 3
        {
            {-1, 0, 0},
            { 1, 0, 0},
            { 0,-1, 0},
            { 0, 1, 0},
            { 0, 0,-1},
            { 0, 0, 1}
        };
        #endif
    };
    
    // Variables (multi-fabs)
    Set::Field<Set::Scalar>    ls_mf;                   // Levelset field variable (current timestep)
    Set::Field<Set::Scalar>    ls_old_mf;               // Levelset field variable (previous timestep)
    Set::Field<Set::Scalar>    XFlux_mf;                // X-direction flux field
    Set::Field<Set::Scalar>    YFlux_mf;                // Y-direction flux field
    Set::Field<Set::Scalar>    ZFlux_mf;                // Z-direction flux field
    Set::Field<Set::Scalar>    Zerols_mf;               // Zero levelset field variable for plotting
    Set::Field<Set::Scalar>    Tube_mf;                 // Narrowband Tube field variable for plotting
    Set::Field<Set::Scalar>    cpt_mf;                  // CPT field variable for plotting
    Set::Field<Set::Scalar>    BA_mf;                   // Narowband Box Array field variable for plotting

    std::unique_ptr<amrex::iMultiFab> Zerols_imf;       // Zero levelset imultifab to track izerols
    std::unique_ptr<amrex::iMultiFab> cpt_imf;          // CPT flag imultifab to track objects
    std::unique_ptr<amrex::iMultiFab> BA_imf;          // CPT flag imultifab to track objects

    // Solver setup and control functions (override)
    void Initialize(int lev) override;
    void TimeStepBegin(Set::Scalar time, int iter) override;
    void TimeStepComplete(Set::Scalar time, int lev) override;
    void Advance(int lev, Set::Scalar time, Set::Scalar dt) override;
    void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;
    void Regrid(int lev, Set::Scalar time) override;

    // Class specific functions
    // Narrowband functions
    void InitializeCPT(int lev, int ls_id);
    void CopyZerolsAndCPTIMFtoMF(int lev);
    void CopyTubeIMFtoMF(int lev, int ls_id);

    void ZeroTheZeros(amrex::iMultiFab& zerols, amrex::iMultiFab& cpt, int ls_id);
    void UpdateNarrowband(int lev, int ls_id);
    //amrex::BoxArray MakeBoxArrayFromIntVects(const amrex::Vector<amrex::IntVect>& ivects, bool simplify = true);
    void ComputeNarrowBandBoxList(int lev, int ls_id);
    void ComputeNarrowBandMapping(int lev, int ls_id);
    void UpdateNarrowBandFlags(int lev, int ls_id);
    void UpdateNarrowbandTubeandMapping(int lev, int ls_id);
    // Geometries
    void ComputeNormal(int lev, int ls_id);
    void ComputeCurvature(int lev, int ls_id);
    void ComputeGeometryQuantities(int lev, int ls_id);
    // Timestepping
    void ComputeAndSetNewTimeStep();
    Set::Scalar GetTimeStep();
    // Advance
    void ClampLS(int lev, int ls_id);
    void Advect(int lev, Set::Scalar time, Set::Scalar dt);
    void UpdateInterfaceVelocity(int lev);
    void Reinitialize(int lev, int ls_id);
    
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
    const int number_of_ghost_cells = 1;     // Number of ghost cells
    
    // Define LevelSetData structure object
    std::vector<LevelSetData> level_sets;
    
    // Define size for narrowband
    const Set::Scalar inner_narrow_band_width = 4.0;
    const Set::Scalar narrow_band_width = 6.0;

    //
    // Definition of user-determined pointer variables.
    //

    IC::IC<Set::Scalar>* ic_ls = nullptr;        // Object used to initialize levelset field
    IC::IC<Set::Scalar>* ic_velocity = nullptr;
    BC::BC<Set::Scalar>* bc_ls = nullptr;        // Object used to update temp levelset field boundary ghost cells
    BC::BC<Set::Scalar>* bc_velocity = nullptr;
    BC::Nothing bc_nothing;
    
    FluxReconstruction reconstruction_method;    // Selected flux reconstruction scheme
    FluxScheme flux_scheme;                      // Selected flux scheme
    TimeSteppingScheme temporal_scheme;          // Selected time-stepping scheme

    Set::Scalar cflNumber = 0.5;                 // CFL number for flow (read from input)
    Set::Scalar cflReinit = 0.4;                 // CFL for reinitialization only

}; //NarrowBandLevelset Class

} // namespace Integrator
#endif
