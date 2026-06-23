//
// Fill a domain based on a list of points to form a closed 2D shape.
// The distance from each point to the nearest edge is found and an error function is then used to compute
// the value of phi between 0 (outside) and 1 (inside) the domain. Only works in 2 dimensions for a single closed polgyon

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

            if (X.back() != X.front())
            // The first and last point must match, if they do not add another point to the end
            {
                X.push_back(X.front());
            }

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Vector x = Set::Position(i, j, k, geom[lev], type);
                Set::Scalar x_cross;
                Set::Scalar min_dist_to_edge = std::numeric_limits<Set::Scalar>::max(); // Mininmum distance from point to nearest edge
                Set::Scalar dist_to_edge; // Distance from point to nearest edge
                Set::Scalar x1;
                Set::Scalar y1;
                Set::Scalar x2;
                Set::Scalar y2;
                int num_cross = 0;

                Set::Scalar px;
                Set::Scalar py;
                Set::Scalar dx;
                Set::Scalar dy;
                Set::Scalar len2;
                Set::Scalar t;
                Set::Vector x_int; // x and y position of the intercept of the mesh point and the edge, assuming edge is infinite

                for (unsigned int n = 0; n < X.size()-1; n++)
                {   

                    x1 = X[n](0);
                    y1 = X[n](1);
                    x2 = X[n+1](0);
                    y2 = X[n+1](1);

                    if ((y1 > x(1)) != (y2 > x(1)))
                    {
                        x_cross = x1 + (x(1) - y1) * (x2 - x1) / (y2 - y1);

                        if (x_cross > x(0))
                            num_cross++;
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
                    px = x(0);
                    py = x(1);
                    dx = x2 - x1;
                    dy = y2 - y1;
                    len2 = dx*dx + dy*dy;
                    t = ((px - x1)*dx + (py - y1)*dy) / len2;

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
        // Diffuseness of the solid boundary
        pp.query_default("eps", value.eps, "0.0", Unit::Length());

        pp.forbid("filename", "use file.name instead");
        // Location of .xy file
        pp.query_file("file.name", filename); 

        // Verbosity (used in parser only)
        pp.query_default("verbose", verbose, 0); 

        // Unitless coordinate multiplier, applied after the rotation
        pp.query_default("mult",value.mult, "1.0", Unit::Less()); 

        // Whether to invert (make IC 1 outside of solid instead of inside solid)
        pp.query_default("invert",value.invert,false); 

        // X-offsed, applied after the rotation
        pp.queryarr_default("x0",value.x0,"0.0 0.0 0.0", Unit::Length());

        // if (pp.contains("x0"))
        //     pp_queryarr("x0",value.x0); // Coordinate offset

        // Center of rigid-body rotation 
        pp.queryarr_default("rotation.center",value.rotation.center,"0.0 0.0 0.0", Unit::Length());
        
        // Amount of rotation about rotation center (clockwise)
        pp.query_default("rotation.angle",value.rotation.angle,"0.0", Unit::Angle());


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

                Set::Scalar cx = value.rotation.center[0];
                Set::Scalar cy = value.rotation.center[1];

                Set::Scalar theta = value.rotation.angle;

                Set::Scalar cos_t = std::cos(theta);
                Set::Scalar sin_t = std::sin(theta);

                Set::Scalar x_shifted = x - cx;
                Set::Scalar y_shifted = y - cy;

                X[0] = cx + x_shifted * cos_t - y_shifted * sin_t;
                X[1] = cy + x_shifted * sin_t + y_shifted * cos_t;

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
    std::vector<Set::Vector> X; // Big X is a vector of vectors of each coordinate read into the file
    Set::Scalar eps;
    Set::Scalar mult = 1.0;
    Set::Vector x0 = Set::Vector::Zero();
    bool invert = false;

    struct{
        Set::Vector center;
        Set::Scalar angle;
    } rotation;
    
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
