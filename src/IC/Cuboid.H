#ifndef IC_CUBOID_H_
#define IC_CUBOID_H_

#include <iostream>
#include <cstdlib>
#include <climits>
#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "mpi.h"
#include "IO/ParmParse.H"

namespace IC
{
class Cuboid : public IC
{
public:
    Cuboid(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    Cuboid(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp.queryclass(name, *this);
    }

    void Define() {
        
    };

    //void Add(const int lev,amrex::Vector<amrex::MultiFab * > &a_field)
    void Add(const int& lev, Set::Field<Set::Scalar>& a_field, Set::Scalar)
    {
        bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

        for (amrex::MFIter mfi(*a_field[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType type = a_field[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");

            amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                Set::Scalar AMREX_D_DECL(x, y, z);
                if (cellcentered)
                {
                    AMREX_D_TERM(x = geom[lev].ProbLo()[0] + ((amrex::Real)(i)+0.5) * geom[lev].CellSize()[0];,
                        y = geom[lev].ProbLo()[1] + ((amrex::Real)(j)+0.5) * geom[lev].CellSize()[1];,
                        z = geom[lev].ProbLo()[2] + ((amrex::Real)(k)+0.5) * geom[lev].CellSize()[2];);
                }
                else
                {
                    AMREX_D_TERM(x = geom[lev].ProbLo()[0] + (amrex::Real)(i)*geom[lev].CellSize()[0];,
                        y = geom[lev].ProbLo()[1] + (amrex::Real)(j)*geom[lev].CellSize()[1];,
                        z = geom[lev].ProbLo()[2] + (amrex::Real)(k)*geom[lev].CellSize()[2];);
                }
                //Set::Scalar rsq = 0.0;


                Set::Scalar xdir = std::abs((x - center[0])) - (length[0] / 2.0);
                Set::Scalar ydir = std::abs((y - center[1])) - (length[1] / 2.0);
            #if AMREX_SPACEDIM==3
                Set::Scalar zdir = std::abs((z - center[2])) - (length[2] / 2.0);
                field(i, j, k, 0) = 0.5 - 0.5 * std::erf((xdir + ydir + zdir) / eps);
                if (field.nComp() > 1) field(i, j, k, 1) = 1. - (0.5 - 0.5 * std::erf((xdir + ydir + zdir) / eps));
            #elif AMREX_SPACEDIM==2
                field(i, j, k, 0) = 0.5 - 0.5 * std::erf((xdir + ydir) / eps);
                if (field.nComp() > 1) field(i, j, k, 1) = 1. - (0.5 - 0.5 * std::erf((xdir + ydir) / eps));
            #endif

                //if (AMREX_D_TERM(std::abs(x - center[0]) <= 0.5 * length[0], && std::abs(y - center[1]) <= 0.5 * length[1], && std::abs(z - center[2]) <= 0.5 * length[2]))
                    //field(i, j, k, 0) = 0.5*(1.0+erf(( x - center[0])/eps)) *0.5*(1.0+erf(( y - center[1])/eps)) ;
                    //"0.5*(1.0+erf(( x - 4.5e-05)/1e-06))* 0.5*(1.0+erf((-x + 5.5e-05)/1e-06))* 0.5*(1.0+erf((y - 4.5e-05)/1e-06))* 0.5*(1.0+erf((-y + 5.5e-05)/1e-06))"
                //else
                   // field(i, j, k, 0) = 0.;
            });
        }

    };

    static void Parse(Cuboid& value, IO::ParmParse& pp)
    {
        pp.queryarr("center", value.center); // Coordinates (X Y Z) of the center of the square/cube. 
        pp.queryarr("length", value.length);
        pp.query("eps", value.eps);
    }

private:
    amrex::Vector<amrex::Real> center;
    amrex::Vector<amrex::Real> length;
    //amrex::Vector<Set::Vector> orientation=(1,0,0);
    Set::Scalar eps;
};
}
#endif
