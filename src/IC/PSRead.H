//
// Fill a domain (region where field=0) with packed spheres (regions where field=1).
// Sphere locations and radii are determined from an xyzr file.
//

#ifndef IC_PSREAD_H_
#define IC_PSREAD_H_

#include "AMReX_Loop.H"
#include "Set/Set.H"
#include "IC/IC.H"
#include "Util/Util.H"
using namespace std;
#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <mpi.h>
#include "IO/ParmParse.H"

namespace IC
{
class PSRead : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "psread";

    enum Type
    {
        Partition,
        Values
    };

    PSRead(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PSRead(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp.queryclass(name, *this);
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
            amrex::LoopConcurrentOnCpu(bx, [=,this] (int i, int j, int k)
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
                if (invert) phi(i,j,k,0) = 1.0 - min_grain_id;
                if (ncomp > 1) phi(i, j, k, 1) = 1.0 - min_grain_id;
            });
        }
    }

    static void Parse(PSRead& value, IO::ParmParse& pp)
    {
        std::string filename;
        int verbose = 0;
        // Diffuseness of the sphere boundary
        pp.query_default("eps", value.eps, "0.0", Unit::Length());

        pp.forbid("filename", "use file.name instead");
        // Location of .xyzr file
        pp.query_file("file.name", filename); 

        // Verbosity (used in parser only)
        pp.query_default("verbose", verbose, 0); 

        // Unitless coordinate multiplier
        pp.query_default("mult",value.mult, "1.0", Unit::Less()); 

        // Whether to invert (1 outside instead of inside spheres)
        pp.query_default("invert",value.invert,false); 

        // X-offset 
        pp.queryarr_default("x0",value.x0,"0 0 0", Unit::Length());

        if (pp.contains("x0"))
            pp_queryarr("x0",value.x0); // Coordinate offset

        Unit unit = Unit::Length();
        // Units of length in the file
        pp.queryunit("file.unit",unit);
        Util::AssertException(  INFO,TEST(unit.isType(Unit::Length())), 
                                "Unit must be of type length but got unit ",unit.normalized_unitstring());
        
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

                Set::Scalar x = (std::stod(strx) * unit).normalized_value();
                Set::Scalar y = (std::stod(stry) * unit).normalized_value();
                #if AMREX_SPACEDIM > 2
                Set::Scalar z = (std::stod(strz) * unit).normalized_value();
                #endif
                Set::Scalar r = (std::stod(strR) * unit).normalized_value();

                Set::Vector X(AMREX_D_DECL(x,y,z));
                X = value.x0 + value.mult*X;

                value.X.push_back(X);
                value.R.push_back(r);
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
    bool invert = false;
};
}
#endif
