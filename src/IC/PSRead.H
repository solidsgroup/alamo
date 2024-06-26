//
// Fill a domain (region where field=0) with packed spheres (regions where field=1).
// Sphere locations and radii are determined from an xyzr file.
//

#ifndef IC_PSREAD_H_
#define IC_PSREAD_H_

#include "Set/Set.H"
#include "IC/IC.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "mpi.h"
#include "IO/ParmParse.H"

namespace IC
{
class PSRead : public IC
{
public:
    enum Type
    {
        Partition,
        Values
    };

    PSRead(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PSRead(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp_queryclass(name, *this);
    }

    void Define() {

    };

    void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar)
    {
        Set::Vector size = Set::Size(geom[lev]);
        int ncomp = a_phi[lev]->nComp();

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType type = a_phi[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                Set::Scalar min_grain_id = 0;

                for (unsigned int n = 0; n < X.size(); n++)
                {
                    Set::Scalar d;
                    Set::Scalar d1 = std::numeric_limits<Set::Scalar>::infinity();

                    d = (x - X[n]).lpNorm<2>();

                    if (geom[0].isPeriodic(0))
                    {
                        d1 = std::min((x - X[n] + size(0) * Set::Vector::Unit(0)).lpNorm<2>(),
                            (x - X[n] - size(0) * Set::Vector::Unit(0)).lpNorm<2>());
                    }
#if AMREX_SPACEDIM > 1
                    if (geom[0].isPeriodic(1))
                    {
                        d1 = std::min((x - X[n] + size(1) * Set::Vector::Unit(1)).lpNorm<2>(),
                            (x - X[n] - size(1) * Set::Vector::Unit(1)).lpNorm<2>());
                    }
#endif
#if AMREX_SPACEDIM > 2
                    if (geom[0].isPeriodic(2))
                    {
                        d1 = std::min((x - X[n] + size(2) * Set::Vector::Unit(2)).lpNorm<2>(),
                            (x - X[n] - size(2) * Set::Vector::Unit(2)).lpNorm<2>());
                    }
#endif

                    d = std::min(d, d1);

                    if (d <= (R[n] + eps))

                    {
                        Set::Scalar m = 0.5 * (1 + erf((-d + R[n]) / (eps)));
                        min_grain_id = min_grain_id + m * (1. - min_grain_id);
                    }
                }

                phi(i, j, k, 0) = min_grain_id;
                if (ncomp > 1) phi(i, j, k, 1) = 1.0 - min_grain_id;
            });
        }
    }

    static void Parse(PSRead& value, IO::ParmParse& pp)
    {
        std::string filename;
        int verbose = 0;
        pp_query("eps", value.eps); // Diffuseness of the sphere boundary
        pp_query_file("filename", filename); // Location of .xyzr file
        pp_query("verbose", verbose); // Verbosity (used in parser only)
        pp_query("mult",value.mult); // Coordinate multiplier
        if (pp.contains("x0"))
            pp_queryarr("x0",value.x0); // Coordinate offset
        std::ifstream datafile(filename);
        std::string line;
        if (datafile.is_open())
        {
            value.X.clear();
            value.R.clear();

            while (getline(datafile, line))
            {
                std::istringstream in(line);

                std::string strx, stry, strz, strR;
                in >> strx >> stry >> strz >> strR;
                value.X.push_back(value.x0 + value.mult*Set::Vector(AMREX_D_DECL(std::stod(strx), std::stod(stry), std::stod(strz))));
                value.R.push_back(value.mult*std::stod(strR));
                if (verbose > 0)
                    Util::Message(INFO, "x=", value.X.back().transpose(), " r=", value.R.back());
            }
            datafile.close();
        }
        else
        {
            Util::Abort(INFO, "Unable to open file ", filename);
        }
    }

private:
    std::vector<Set::Vector> X;
    std::vector<Set::Scalar> R;
    Set::Scalar eps;
    Set::Scalar mult = 1.0;
    Set::Vector x0 = Set::Vector::Zero();
};
}
#endif
