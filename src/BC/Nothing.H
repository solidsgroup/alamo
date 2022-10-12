#ifndef BC_NOTHING_H_
#define BC_NOTHING_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>

namespace BC
{
class Nothing
    : public BC<Set::Scalar>
{

public:
    Nothing () {};
    virtual void FillBoundary (amrex::BaseFab<Set::Scalar>&,
                    const amrex::Box &,
                    int, int, int,
                    amrex::Real,
                    Orientation /*face*/ = Orientation::All,
                    const amrex::Mask * /*mask*/ = nullptr) override {};

    virtual amrex::BCRec GetBCRec() override {return amrex::BCRec();}
};
}
#endif
