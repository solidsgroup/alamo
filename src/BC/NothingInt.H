#ifndef BC_NOTHING_INT_H_
#define BC_NOTHING_INT_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>


#include "Set/Set.H"
#include "BC/BC.H"

namespace BC {

class NothingInt : public BC<Set::IntScalar> {
public:
    NothingInt() = default;
    virtual ~NothingInt() {}

    /// Override FillBoundary for integer fields
    void FillBoundary(amrex::BaseFab<Set::IntScalar>& in,
                      const amrex::Box& box,
                      int ngrow, int dcomp, int ncomp, amrex::Real time,
                      Orientation face = Orientation::All,
                      const amrex::Mask* mask = nullptr) override;

    /// Override GetBCRec
    amrex::BCRec GetBCRec() override {
        return amrex::BCRec();
    }
};


} // namespace BC

#endif
