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

        PSRead(amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}

        PSRead(amrex::Vector<amrex::Geometry> &_geom, IO::ParmParse &pp, std::string name) : IC(_geom) 
        {
            pp.queryclass(name,*this);
        }

        void Define(){

        };

        void Add(const int &lev, Set::Field<Set::Scalar> &a_phi)
        {
            Set::Vector size;
            AMREX_D_TERM(size(0) = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];,
                         size(1) = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];,
                         size(2) = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2];)

            for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.tilebox();
                bx.grow(a_phi[lev]->nGrow());
                amrex::Array4<Set::Scalar> const &phi = a_phi[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                   {
                                       Set::Vector x;
                                       Set::Scalar min_grain_id = 0;

                                       AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
                                                    x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
                                                    x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);

                                       for (unsigned int n = 0; n < X.size(); n++)
                                       {
                                           Set::Scalar d;
                                           Set::Scalar d1;

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

                                       phi(i, j, k) = min_grain_id;
                                   });
            }
        }

        static void Parse(PSRead &value, IO::ParmParse &pp)
        {
            std::string filename;
            pp.query("eps",value.eps);
            pp.query("filename", filename);
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
                    value.X.push_back(Set::Vector(AMREX_D_DECL(std::stod(strx),std::stod(stry),std::stod(strz))));
                    value.R.push_back(std::stod(strR));
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
    };
}
#endif
