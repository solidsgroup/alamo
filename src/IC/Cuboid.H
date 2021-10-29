#ifndef IC_CUBOID_H_
#define IC_CUBOID_H_

#include <iostream>
#include <cstdlib>
#include <climits>
#include "IC/IC.H"
#include "Util/Util.H"

/// \class Sphere
/// \brief Initialize a spherical inclusion
namespace IC
{
class Cuboid : public IC
{
public:
    Cuboid (amrex::Vector<amrex::Geometry> &_geom) :
        IC(_geom)
    {
        amrex::ParmParse pp("ic");
        pp.queryarr("center",center);
        pp.queryarr("length",length);
    }
  
    void Add(const int lev,
            amrex::Vector<amrex::MultiFab * > &a_field)
    {
        bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

        for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            bx.grow(a_field[lev]->nGrow());

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

                    Set::Scalar AMREX_D_DECL(x,y,z);
                    if (cellcentered)
                    {
                        AMREX_D_TERM(x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                                    y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                                    z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                    }
                    else
                    {
                        AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i) * geom[lev].CellSize()[0];,
                                    y = geom[lev].ProbLo()[1] + (amrex::Real)(j) * geom[lev].CellSize()[1];,
                                    z = geom[lev].ProbLo()[2] + (amrex::Real)(k) * geom[lev].CellSize()[2];);
                    }
                    //Set::Scalar rsq = 0.0;

                    if (AMREX_D_TERM(std::abs(x-center[0]) <= 0.5*length[0], && std::abs(y-center[1]) <= 0.5*length[1], && std::abs(z-center[2]) <= 0.5*length[2]))
                        field(i,j,k,0) = 0.;
                    else
                        field(i,j,k,0) = 1.;  
                });
        }

    };
  
private:
    amrex::Vector<amrex::Real> center;
    amrex::Vector<amrex::Real> length;
};
}
#endif
