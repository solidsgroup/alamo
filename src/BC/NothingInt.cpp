#include "BC/NothingInt.H"
#include "IO/ParmParse.H" // Including IO for ParmParse

namespace BC {

/// \brief FillBoundary implementation for NothingInt
void NothingInt::FillBoundary(amrex::BaseFab<Set::IntScalar>& in,
                              const amrex::Box& box,
                              int ngrow, int dcomp, int ncomp, amrex::Real time,
                              Orientation face, const amrex::Mask* mask) {
    // No action is taken; this is a "do-nothing" boundary condition.
}




} // namespace BC
