//
// Fill a domain based on a list of points to form a closed 2D shape.
// The distance from each point to the nearest edge is found and an error function is then used to compute
// the vlude of phi between 0 (outside) and 1 (inside).
//

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
            double small = 1e-8;
            m_edge.resize(X.size() - 1);
            b_edge.resize(X.size() - 1);

            for (unsigned int n = 0; n < X.size()-1; n++)
                {   
                    // Compute the slope and y intercept of each edge
                    m_edge[n] = (X[n + 1](1) - X[n](1))/(X[n + 1](0) - X[n](0) + small);
                    b_edge[n] = (X[n](1) - m_edge[n]*X[n](0));
                }

            // std::cout << "X.size() = " << X.size() << std::endl;
            // std::cout << "m_edge.size() = " << m_edge.size() << std::endl;
            // std::cout << "b_edge.size() = " << b_edge.size() << std::endl;

            // m_edge.resize(X.size() - 1);
            // b_edge.resize(X.size() - 1);

            // for (unsigned int n = 0; n < X.size() - 1; n++)
            // {
            //     std::cout << "\nIteration n = " << n << std::endl;

            //     std::cout << "Accessing X[" << n << "] and X[" << n + 1 << "]" << std::endl;

            //     std::cout << "X[" << n << "] = ("
            //             << X[n](0) << ", "
            //             << X[n](1) << ")" << std::endl;

            //     std::cout << "X[" << n + 1 << "] = ("
            //             << X[n + 1](0) << ", "
            //             << X[n + 1](1) << ")" << std::endl;

            //     double dx = X[n + 1](0) - X[n](0) + small;
            //     double dy = X[n + 1](1) - X[n](1);

            //     std::cout << "dx = " << dx << ", dy = " << dy << std::endl;

            //     m_edge[n] = dy / dx;
            //     b_edge[n] = X[n](1) - m_edge[n] * X[n](0);

            //     std::cout << "m_edge[" << n << "] = " << m_edge[n] << std::endl;
            //     std::cout << "b_edge[" << n << "] = " << b_edge[n] << std::endl;
            // }

            // std::cout << "Finished computing slopes \n";

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
                    Set::Scalar m_r; // Slope of the line normal to the edge
                    Set::Scalar b_r; // Y-intercept of the line normal to the edge
                    Set::Vector x_int; // x and y position of the intercept of the mesh point and the edge, assuming edge is infinite
                    // std::cout << "X[n](1) = " << X[n](1) << ",X[n][1] = " << X[n][1] << std::endl;
                    if ((X[n](1) > x(1) && X[n+1](1) > x(1)) || (X[n](1) < x(1) && X[n+1](1) < x(1)))
                    {
                        num_cross += 0;
                    } else
                    {
                        // x_cross = (x(1) - b_edge[n])/(m_edge[n]);
                        // std::cout << "\nComputing x_cross for edge n = " << n << std::endl;
                        // std::cout << "x(1) = " << x(1) << std::endl;
                        // std::cout << "b_edge[" << n << "] = " << b_edge[n] << std::endl;
                        // std::cout << "m_edge[" << n << "] = " << m_edge[n] << std::endl;

                        double numerator = x(1) - b_edge[n];

                        // std::cout << "numerator = x(1) - b_edge[" << n << "] = " << numerator << std::endl;

                        if (std::abs(m_edge[n]) < 1e-12)
                        {
                            std::cout << "WARNING: m_edge[" << n
                                    << "] is near zero!" << std::endl;
                        }

                        x_cross = numerator / m_edge[n];

                        // std::cout << "x_cross = " << x_cross << std::endl;

                        if (x_cross > x[0])
                        {   // If the horxiontal ray crosses the intercept, add 1 to n_cross
                            num_cross += 1;
                        } else
                        {
                            num_cross += 0;
                        }
                    }

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

                    // m_r = -1/(m_edge[n]+small);
                    // b_r = x(1) - m_r*x(0);

                    // x_int(0) = (b_edge[n] - b_r)/(m_r - m_edge[n] + small);
                    // x_int(1) = m_edge[n]*x_int(0) + b_edge[n];

                    // if ((x_int(0) > X[n](0)) && (x_int(0) < X[n+1](0))) {

                    // }

                    // if (((x_int(0) - X[n](0)) > 0 && (x_int(0) - X[n+1](0) < 0)) && ((x_int(1) - X[n](1)) > 0 && (x_int(1) - X[n+1](1) < 0)))
                    // {
                    //     dist_to_edge = dist(x(0), x(1), x_int(0), x_int(1));
                    // } else if (((x_int(0) - X[n](0)) < 0 && (x_int(0) - X[n+1](0) > 0)) && ((x_int(1) - X[n](1)) < 0 && (x_int(1) - X[n+1](1) > 0)))
                    // {
                    //     dist_to_edge = dist(x(0), x(1), x_int(0), x_int(1));
                    // } else
                    // {
                    //     dist_to_edge = std::min(dist(x(0), x(1), X[n](0), X[n](1)), dist(x(0), x(1), X[n+1](0), X[n+1](1)));
                    // }
                    min_dist_to_edge = std::min(min_dist_to_edge, dist_to_edge);

                }
                if (num_cross % 2 == 0) 
                {
                    min_dist_to_edge = min_dist_to_edge; // If the number of crossings is even, the point is inside the shape
                } else 
                {
                    min_dist_to_edge = -1*min_dist_to_edge;
                }

                double denom = std::sqrt(2.0) * eps;
                double arg   = min_dist_to_edge / denom;
                double tanh_val = std::tanh(arg);
                double phi_val  = 0.5 * (1.0 - tanh_val);

                // std::cout << "i,j,k = " << i << ", " << j << ", " << k << std::endl;

                // std::cout << "num_cross = "<< num_cross << std::endl;

                // std::cout << "min_dist_to_edge = " << min_dist_to_edge << std::endl;

                // std::cout << "phi_val = " << phi_val << "\n" << std::endl;
                
                phi(i,j,k) = 0.5 * (1.0 - std::tanh(min_dist_to_edge /(std::sqrt(2.0) * eps)));
                
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
