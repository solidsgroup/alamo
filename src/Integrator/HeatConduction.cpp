#include "HeatConduction.H"
#include "BC/Constant.H"
#include "IC/Constant.H"

namespace Integrator
{
/// \fn    HeatConduction::Integrator::Integrator
///
/// Read in the following simulation parameters
///
///     heat.alpha                (default 1.0)
///     heat.refinement_threshold (default 0.01)
///     ic.type
///

/// Initialize initial condition pointer #ic, and register
/// the #Temp, #Temp_old Multifab arrays.

//HeatConduction::HeatConduction() :
//    Integrator()

/// \fn HeatConduction::Integrator::Initialize
///
/// Use the #ic object to initialize #Temp
//void
//HeatConduction::Initialize (int lev)


/// \fn    Integrator::HeatConduction::Advance
///
/// Integrate the heat diffusion equation
/// \f[\nabla^2T = \alpha \frac{\partial T}{\partial t}\f]
/// using an explicit forward Euler method.
/// \f$\alpha\f$ is stored in #alpha
//void
//HeatConduction::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)


/// \fn    Integrator::HeatConduction::TagCellsForRefinement
///
/// The following criterion is used to determine if a cell should be refined:
/// \f[|\nabla T|\,|\mathbf{r}| > h\f]
/// where
/// \f[\mathbf{r} = \sqrt{\Delta x_1^2 + \Delta x_2^2 + \Delta x_3^2}\f]
/// and \f$h\f$ is stored in #refinement_threshold
//void
//HeatConduction::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)

}
