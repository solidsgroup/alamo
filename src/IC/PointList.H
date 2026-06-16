//
// Fill a domain based on a list of points to form a closed 2D shape.
// The distance from each point to the nearest edge is found and an error function is then used to compute
// the vlude of phi between 0 (outside) and 1 (inside). Only works in 2 dimensions

#ifndef IC_POINTLIST_H_
#define IC_POINTLIST_H_

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
class PointList : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "pointlist";

    enum Type
    {
        Partition,
        Values
    };

    PointList(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    PointList(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name) : IC(_geom)
    {
        pp.queryclass(name, *this);
    }
    void Define() {

    };

    void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar)
    {   
        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {   
            if (AMREX_SPACEDIM != 2)
            {
                amrex::Abort("This code only supports 2D (AMREX_SPACEDIM must be 2)");
            }

            amrex::Box bx;
            amrex::IndexType type = a_phi[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
            else Util::Abort(INFO, "Unkonwn index type");
            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);
            double small = 1e-8;
            m_edge.resize(X.size() - 1);
            b_edge.resize(X.size() - 1);

            if (X.back() != X.front())
            // The first and last point must match, if they do not add another point to the end
            {
                X.push_back(X.front());
            }
            for (unsigned int n = 0; n < X.size()-1; n++)
                {   
                    // Compute the slope and y intercept of each edge
                    m_edge[n] = (X[n + 1](1) - X[n](1))/(X[n + 1](0) - X[n](0) + small);
                    b_edge[n] = (X[n](1) - m_edge[n]*X[n](0));
                }

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                Set::Scalar x_cross;
                Set::Scalar min_dist_to_edge = std::numeric_limits<Set::Scalar>::max(); // Mininmum distance from point to nearest edge
                Set::Scalar dist_to_edge; // Distance from point to nearest edge
                int num_cross = 0;

                // std::cout << "X.size() = " << X.size() << std::endl;
                // std::cout << "X[0] = " << X[0] << "X[1] = " << X[1] << "X[2] = " << X[2] << std::endl;

                for (unsigned int n = 0; n < X.size()-1; n++)
                {   
                    Set::Vector x_int; // x and y position of the intercept of the mesh point and the edge, assuming edge is infinite
                    // std::cout << "X[n](1) = " << X[n](1) << ",X[n][1] = " << X[n][1] << std::endl;
                    if ((X[n](1) > x(1) && X[n+1](1) > x(1)) || (X[n](1) < x(1) && X[n+1](1) < x(1)))
                    {
                        num_cross += 0;
                    } else
                    {
                        double numerator = x(1) - b_edge[n];

                        x_cross = numerator / m_edge[n];

                        if (x_cross > x[0])
                        {   // If the horxiontal ray crosses the intercept, add 1 to n_cross
                            num_cross += 1;
                        } else
                        {
                            num_cross += 0;
                        }
                    }

                    // Computes the closest point on a line segment (X[n] → X[n+1]) to a query point x,
                    // using vector projection onto the segment.
                    //
                    // Let A = X[n], B = X[n+1], and P = x.
                    //
                    // The parameter t is the normalized projection of P onto the infinite line AB:
                    //
                    //     t = ((P - A) · (B - A)) / |B - A|^2
                    //
                    // If t ∈ [0, 1], the perpendicular projection lies within the segment.
                    // If t < 0, the closest point is A.
                    // If t > 1, the closest point is B.
                    //
                    // The resulting closest point x_int is then used to compute the Euclidean
                    // distance from P to the segment, and the minimum distance over all edges is tracked.

                    double x1 = X[n](0);
                    double y1 = X[n](1);
                    double x2 = X[n+1](0);
                    double y2 = X[n+1](1);
                    double px = x(0);
                    double py = x(1);
                    double dx = x2 - x1;
                    double dy = y2 - y1;
                    double len2 = dx*dx + dy*dy;
                    double t = ((px - x1)*dx + (py - y1)*dy) / len2;

                    if (t < 0.0)
                    {
                        x_int(0) = x1;
                        x_int(1) = y1;
                    }
                    else if (t > 1.0)
                    {
                        x_int(0) = x2;
                        x_int(1) = y2;
                    }
                    else
                    {
                        x_int(0) = x1 + t*dx;
                        x_int(1) = y1 + t*dy;
                    }

                    dist_to_edge = dist(x(0), x(1), x_int(0), x_int(1));
                    min_dist_to_edge = std::min(min_dist_to_edge, dist_to_edge);

                }
                if (num_cross % 2 == 0) 
                {
                    min_dist_to_edge = 1*min_dist_to_edge; // If the number of crossings is even, the point is inside the shape
                } else 
                {
                    min_dist_to_edge = -1*min_dist_to_edge;
                }
                
                phi(i,j,k) = 0.5 * (1.0 - std::tanh(min_dist_to_edge /(std::sqrt(2.0) * eps)));

                if (invert) {
                    phi(i,j,k) = 1 - phi(i,j,k);
                }
                
            });
        }
    }

    static void Parse(PointList& value, IO::ParmParse& pp)
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

                Set::Vector X(AMREX_D_DECL(x,y,z));
                X = value.x0 + value.mult*X;

                value.X.push_back(X);
                if (verbose > 0)
                    Util::Message(INFO, "x=", value.X.back().transpose());
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
    std::vector<Set::Scalar> m_edge;
    std::vector<Set::Scalar> b_edge;
    bool invert = false;
    
    double dist(double x1, double y1, double x2, double y2)
    {   
        // Find distance between 2 points
        double dx = x2 - x1;
        double dy = y2 - y1;

        return std::sqrt(dx*dx + dy*dy);
    }
};
}
#endif
