//
// Fill a domain based on a list of points to form one or more closed 2D shapes.
// The distance from each point to the nearest edge is found and an error function is then used to compute
// the value of phi between 0 (outside) and 1 (inside) the domain. Only works in 2 dimensions.
// Points belonging to different polygons are distinguished using the ObjNum column of the input file:
// whenever ObjNum changes by more than obj_num_threshold between consecutive lines, a new polygon begins.
// Each polygon's field is computed into phi_temp. Polygons are then combined into either the "solid" field
// (pointwise max over all non-void polygons) or the "void" field (pointwise max over all void polygons,
// as selected by the per-polygon invert list). The void field is inverted (zero inside, one outside the
// void polygons) and multiplied into the solid field, punching holes wherever a void polygon overlaps it.
// invert_all is then applied once at the end to flip the whole finished field if desired.

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
        if (AMREX_SPACEDIM != 2)
        {
            amrex::Abort("This code only supports 2D (AMREX_SPACEDIM must be 2)");
        }

        if ((int)phi_temp.size() != a_phi.size()) phi_temp.resize(a_phi.size());
        if ((int)phi_void.size() != a_phi.size()) phi_void.resize(a_phi.size());
        if (!phi_temp[lev] ||
            phi_temp[lev]->boxArray() != a_phi[lev]->boxArray() ||
            phi_temp[lev]->DistributionMap() != a_phi[lev]->DistributionMap())
        {
            phi_temp[lev] = std::make_unique<amrex::MultiFab>(a_phi[lev]->boxArray(), a_phi[lev]->DistributionMap(),
                                                                a_phi[lev]->nComp(), a_phi[lev]->nGrow());
        }
        if (!phi_void[lev] ||
            phi_void[lev]->boxArray() != a_phi[lev]->boxArray() ||
            phi_void[lev]->DistributionMap() != a_phi[lev]->DistributionMap())
        {
            phi_void[lev] = std::make_unique<amrex::MultiFab>(a_phi[lev]->boxArray(), a_phi[lev]->DistributionMap(),
                                                                a_phi[lev]->nComp(), a_phi[lev]->nGrow());
        }

        bool first_solid = true;
        bool first_void = true;
        bool any_void = false;

        for (unsigned int p = 0; p < polygons.size(); p++)
        {
            std::vector<Set::Vector> polygon = polygons[p];
            if (polygon.back() != polygon.front())
            // The first and last point must match, if they do not add another point to the end
            {
                polygon.push_back(polygon.front());
            }
            bool is_void = (p < invert.size() && invert[p] != 0);
            bool is_first = is_void ? first_void : first_solid;

            for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx;
                amrex::IndexType type = a_phi[lev]->ixType();
                if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
                else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
                else Util::Abort(INFO, "Unkonwn index type");
                amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const& phi_t = phi_temp[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const& phi_v = phi_void[lev]->array(mfi);

                std::vector<Set::Vector> X = polygon;

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

                    phi_t(i,j,k) = 0.5 * (1.0 - std::tanh(min_dist_to_edge /(std::sqrt(2.0) * eps)));

                    // Combine this polygon's field into either the running solid union or the running
                    // void union (pointwise max), depending on whether it is marked as a void polygon.
                    if (is_void) phi_v(i,j,k) = is_first ? phi_t(i,j,k) : std::max(phi_v(i,j,k), phi_t(i,j,k));
                    else         phi(i,j,k)   = is_first ? phi_t(i,j,k) : std::max(phi(i,j,k), phi_t(i,j,k));
                });
            }

            if (is_void) { first_void = false; any_void = true; }
            else         first_solid = false;
        }

        if (any_void)
        // Punch holes into the solid field wherever a void polygon overlaps it: phi_void is zero inside
        // the void polygons and one outside them, so multiplying it in zeros out the void regions.
        {
            for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx;
                amrex::IndexType type = a_phi[lev]->ixType();
                if (type == amrex::IndexType::TheCellType())      bx = mfi.growntilebox();
                else if (type == amrex::IndexType::TheNodeType()) bx = mfi.grownnodaltilebox();
                else Util::Abort(INFO, "Unkonwn index type");
                amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const& phi_v = phi_void[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    phi(i,j,k) = phi(i,j,k) * (1.0 - phi_v(i,j,k));
                });
            }
        }

        if (invert_all)
        {
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
                    phi(i,j,k) = 1 - phi(i,j,k);
                });
            }
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

        // Whether to invert the finished field as a whole (make IC 1 outside of solid instead of inside solid)
        pp.query_default("invert_all",value.invert_all,false);

        // X-offsed, applied after the rotation
        pp.queryarr_default("x0",value.x0,"0.0 0.0 0.0", Unit::Length());

        // Center of rigid-body rotation 
        pp.queryarr_default("rotation.center",value.rotation.center,"0.0 0.0 0.0", Unit::Length());
        
        // Amount of rotation about rotation center (clockwise)
        pp.query_default("rotation.angle",value.rotation.angle,"0.0", Unit::Angle());

        // ObjNum jump (between consecutive points in the file) larger than this starts a new polygon
        pp.query_default("obj_num_threshold",value.objnum_threshold,"0.5", Unit::Less());


        Unit unit = Unit::Length();
        // Units of length in the file
        pp.queryunit("file.unit",unit);
        Util::AssertException(  INFO,TEST(unit.isType(Unit::Length())), 
                                "Unit must be of type length but got unit ",unit.normalized_unitstring());
        
        std::ifstream datafile(filename);
        std::string line;
        if (datafile.is_open())
        {
            value.polygons.clear();
            bool have_prev_obj_num = false;
            Set::Scalar prev_obj_num = 0.0;

            while (getline(datafile, line))
            {
                std::istringstream in(line);

                std::string strx, stry, strz, strObjNum;
                in >> strx >> stry >> strz;

                Set::Scalar x = (std::stod(strx) * unit).normalized_value();
                Set::Scalar y = (std::stod(stry) * unit).normalized_value();
                #if AMREX_SPACEDIM > 2
                Set::Scalar z = (std::stod(strz) * unit).normalized_value();
                #endif

                // ObjNum is optional; if the file has no fourth column, treat every point as one polygon
                Set::Scalar ObjNum = have_prev_obj_num ? prev_obj_num : 0.0;
                if (in >> strObjNum) ObjNum = std::stod(strObjNum);

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

                // A jump in ObjNum larger than the threshold marks the start of a new polygon
                if (!have_prev_obj_num || std::abs(ObjNum - prev_obj_num) > value.objnum_threshold)
                {
                    value.polygons.push_back(std::vector<Set::Vector>());
                }
                have_prev_obj_num = true;
                prev_obj_num = ObjNum;

                value.polygons.back().push_back(X);
                if (verbose > 0)
                    Util::Message(INFO, "x=", value.polygons.back().back().transpose());
            }
            datafile.close();
        }
        else
        {
            Util::Abort(INFO, "Unable to open file ", filename);
        }

        // Per-polygon void flag: nonzero marks a polygon as a void/negative region that is cut out of
        // the solid field instead of being unioned into it. Must have one entry per polygon if given.
        if (pp.queryarr("invert", value.invert))
        {
            Util::AssertException(INFO, TEST(value.invert.size() == value.polygons.size()),
                "'invert' must have one entry per polygon: expected ", value.polygons.size(),
                " but got ", value.invert.size());
        }
        else
        {
            value.invert.assign(value.polygons.size(), 0);
        }
    }

private:
    std::vector<std::vector<Set::Vector>> polygons; // Each entry is the ordered list of vertices for one closed polygon
    Set::Scalar eps;
    Set::Scalar mult = 1.0;
    Set::Vector x0 = Set::Vector::Zero();
    Set::Field<Set::Scalar> phi_temp;
    Set::Field<Set::Scalar> phi_void;
    Set::Scalar objnum_threshold = 0.5;
    std::vector<int> invert; // per-polygon void flag: nonzero = void/negative region
    bool invert_all = false;

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
