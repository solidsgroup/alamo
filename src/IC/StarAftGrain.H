//
// Initial condition for a star-aft grain geometry, based on NAWC Motor No. 6.
// The grain consists of a cylindrical (tube) bore at the head end transitioning
// to a star-shaped port at the aft end. Cross-section parameters are specified
// at discrete axial stations and linearly interpolated between them.
//
// Geometry is defined in the (x1=x, x2=y, x3=z) coordinate system where
// x3 is the axial direction. The motor axis lies along x3.
//
// Station parameters (from Willcox et al. 2007, Fig. 7, NAWC Motor No. 6):
//   W    = web thickness
//   R1   = inner port radius (bore radius)
//   R3   = star tip radius (outer radius of star points)
//   R4   = fillet radius at star point tips
//   nfin = number of star fins (points)
//   fin_half_angle = half-angle of each fin (degrees)
//
// phi = 1 inside propellant, phi = 0 inside gas port.
// A smooth tanh transition of width eps is applied at the interface.
//

#ifndef IC_STARAFTGRAIN_H_
#define IC_STARAFTGRAIN_H_

#include <cmath>
#include <vector>
#include <algorithm>

#include "IO/ParmParse.H"
#include "IC/IC.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace IC
{

class StarAftGrain : public IC<Set::Scalar>
{
public:
    static constexpr const char* name = "staraftgrain";

    StarAftGrain(amrex::Vector<amrex::Geometry>& _geom) : IC(_geom) {}

    StarAftGrain(amrex::Vector<amrex::Geometry>& _geom, IO::ParmParse& pp, std::string name)
        : IC(_geom)
    {
        pp.queryclass(name, *this);
    }

    // -----------------------------------------------------------------------
    // Per-station cross-section parameters
    // -----------------------------------------------------------------------
    struct Station
    {
        Set::Scalar z;              // Axial position [m]
        Set::Scalar R1;             // Inner bore / port radius [m]
        Set::Scalar R3;             // Star tip outer radius [m]  (=R1 for tube stations)
        Set::Scalar R4;             // Fillet radius at fin tips [m] (0 for tube stations)
        bool        is_star;        // True if this station has star cross-section
    };

    // -----------------------------------------------------------------------
    // Signed distance function for a star cross-section at a given (x,y)
    //
    // Returns positive values inside the propellant web,
    //         negative values inside the gas port.
    //
    // The star is defined by nfin fins equally spaced in angle.
    // Each fin is a wedge of half-angle `half_angle` (radians), with:
    //   R1  = root radius (inner bore circle)
    //   R3  = tip radius  (outer tip circle)
    //   R4  = tip fillet radius
    // -----------------------------------------------------------------------
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar starSDF(
        Set::Scalar x, Set::Scalar y,
        Set::Scalar R1, Set::Scalar R3, Set::Scalar R4,
        int nfin, Set::Scalar half_angle)
    {
        Set::Scalar r     = std::sqrt(x * x + y * y);
        Set::Scalar theta = std::atan2(y, x);

        // Angular period of one fin
        Set::Scalar period = 2.0 * Set::Constant::Pi / (Set::Scalar)nfin;

        // Map theta into [0, period) and find angular distance to nearest fin center
        Set::Scalar theta_mod = std::fmod(theta + 4.0 * Set::Constant::Pi, period);
        Set::Scalar dtheta    = std::min(theta_mod, period - theta_mod);

        // Effective port radius at this angle:
        //   - Within the fin half-angle: port opens out to R3 (star tip)
        //   - Outside the fin: port is just the bore circle R1
        // We interpolate smoothly between R1 and R3 based on angular position.
        Set::Scalar t         = std::min(dtheta / half_angle, 1.0);
        Set::Scalar R_port    = R1 + (R3 - R4 - R1) * (1.0 - t * t);

        // SDF: positive = propellant (outside port), negative = gas (inside port)
        return r - R_port;
    }

    // -----------------------------------------------------------------------
    // Simple tube SDF: positive outside bore (propellant), negative inside
    // -----------------------------------------------------------------------
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar tubeSDF(Set::Scalar x, Set::Scalar y, Set::Scalar R1)
    {
        return std::sqrt(x * x + y * y) - R1;
    }

    // -----------------------------------------------------------------------
    // Main IC fill loop
    // -----------------------------------------------------------------------
    void Add(const int& lev, Set::Field<Set::Scalar>& a_phi, Set::Scalar) override
    {
        Util::Assert(INFO, TEST(stations.size() >= 2),
                     "StarAftGrain: need at least 2 stations");

        // Copy station data into flat arrays for GPU lambda capture
        int nstations = (int)stations.size();
        std::vector<Set::Scalar> st_z(nstations), st_R1(nstations);
        std::vector<Set::Scalar> st_R3(nstations), st_R4(nstations);
        std::vector<int>         st_star(nstations);

        for (int s = 0; s < nstations; s++)
        {
            st_z[s]    = stations[s].z;
            st_R1[s]   = stations[s].R1;
            st_R3[s]   = stations[s].R3;
            st_R4[s]   = stations[s].R4;
            st_star[s] = stations[s].is_star ? 1 : 0;
        }

        // AMReX GPU-accessible copies
        // Use simple host arrays - compatible with CPU builds
        const Set::Scalar* p_z    = st_z.data();
        const Set::Scalar* p_R1   = st_R1.data();
        const Set::Scalar* p_R3   = st_R3.data();
        const Set::Scalar* p_R4   = st_R4.data();
        const int*         p_star = st_star.data();

        const Set::Scalar  _eps        = eps;
        const Set::Scalar  _half_angle = half_angle;
        const int          _nfin       = nfin;
        const int          _nstations  = nstations;
        const Set::Scalar  _z_motor_lo = z_motor_lo;
        const Set::Scalar  _R_case     = R_case;

        for (amrex::MFIter mfi(*a_phi[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx;
            amrex::IndexType type = a_phi[lev]->ixType();
            if (type == amrex::IndexType::TheCellType())
                bx = mfi.growntilebox();
            else if (type == amrex::IndexType::TheNodeType())
                bx = mfi.grownnodaltilebox();
            else
                Util::Abort(INFO, "StarAftGrain: unknown index type");

            amrex::Array4<Set::Scalar> const& phi = a_phi[lev]->array(mfi);


            Util::Message(INFO, "StarAftGrain: entering Add() loop, lev=", lev);
            Util::Message(INFO, "StarAftGrain: exiting Add() loop");
            
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                // Physical coordinates: x1,x2 are transverse; x3 is axial
                Set::Scalar x1, x2, x3;
                if (type == amrex::IndexType::TheCellType())
                {
                    x1 = geom[lev].ProbLo()[0] + ((Set::Scalar)i + 0.5) * geom[lev].CellSize()[0];
                    x2 = geom[lev].ProbLo()[1] + ((Set::Scalar)j + 0.5) * geom[lev].CellSize()[1];
#if AMREX_SPACEDIM > 2
                    x3 = geom[lev].ProbLo()[2] + ((Set::Scalar)k + 0.5) * geom[lev].CellSize()[2];
#else
                    x3 = _z_motor_lo; // 2D fallback: use head-end station
#endif
                }
                else
                {
                    x1 = geom[lev].ProbLo()[0] + (Set::Scalar)i * geom[lev].CellSize()[0];
                    x2 = geom[lev].ProbLo()[1] + (Set::Scalar)j * geom[lev].CellSize()[1];
#if AMREX_SPACEDIM > 2
                    x3 = geom[lev].ProbLo()[2] + (Set::Scalar)k * geom[lev].CellSize()[2];
#else
                    x3 = _z_motor_lo;
#endif
                }

                // Transverse position relative to motor axis
                Set::Scalar x = x1;
                Set::Scalar y = x2;

                // -----------------------------------------------------------
                // Find the two bracketing stations for axial position x3
                // -----------------------------------------------------------
                int s_lo = 0;
                int s_hi = _nstations - 1;

                // Clamp to motor extent
                if (x3 <= p_z[0])
                {
                    s_lo = 0; s_hi = 0;
                }
                else if (x3 >= p_z[_nstations - 1])
                {
                    s_lo = _nstations - 1; s_hi = _nstations - 1;
                }
                else
                {
                    for (int s = 0; s < _nstations - 1; s++)
                    {
                        if (x3 >= p_z[s] && x3 < p_z[s + 1])
                        {
                            s_lo = s;
                            s_hi = s + 1;
                            break;
                        }
                    }
                }

                // Interpolation weight
                Set::Scalar t_interp = 0.0;
                if (s_hi != s_lo)
                    t_interp = (x3 - p_z[s_lo]) / (p_z[s_hi] - p_z[s_lo]);

                // Linearly interpolated cross-section parameters
                Set::Scalar R1 = p_R1[s_lo] + t_interp * (p_R1[s_hi] - p_R1[s_lo]);
                Set::Scalar R3 = p_R3[s_lo] + t_interp * (p_R3[s_hi] - p_R3[s_lo]);
                Set::Scalar R4 = p_R4[s_lo] + t_interp * (p_R4[s_hi] - p_R4[s_lo]);
                int         is_star = (p_star[s_lo] || p_star[s_hi]) ? 1 : 0;

                // -----------------------------------------------------------
                // Compute signed distance to the port boundary
                // Positive SDF = propellant, Negative SDF = gas port
                // -----------------------------------------------------------
                Set::Scalar sdf_port;
                if (is_star)
                    sdf_port = starSDF(x, y, R1, R3, R4, _nfin, _half_angle);
                else
                    sdf_port = tubeSDF(x, y, R1);

                // Also enforce the outer case boundary (propellant doesn't exist outside case)
                Set::Scalar sdf_case = _R_case - std::sqrt(x * x + y * y);

                // Take the minimum: inside propellant only if inside case AND outside port
                Set::Scalar sdf = std::min(sdf_port, sdf_case);

                // Smooth tanh transition
                Set::Scalar eta;
                if (_eps > 0.0)
                    eta = 0.5 * (1.0 + std::tanh(sdf / _eps));
                else
                    eta = (sdf >= 0.0) ? 1.0 : 0.0;

                phi(i, j, k, 0) = eta;
            });
        }
    }

    // -----------------------------------------------------------------------
    // Input file parser
    // -----------------------------------------------------------------------
    static void Parse(StarAftGrain& value, IO::ParmParse& pp)
    {
        // Interface diffuseness [length units]
        pp.query_default("eps", value.eps, "0.001", Unit::Length());

        // Outer case radius [length units]
        pp.query_default("R_case", value.R_case, "0.127", Unit::Length());

        // Motor head-end axial coordinate (used as 2D fallback)
        pp.query_default("z_motor_lo", value.z_motor_lo, "0.0", Unit::Length());

        // Number of star fins
        pp.query_default("nfin", value.nfin, 6);

        // Fin half-angle [degrees] — converted to radians
        Set::Scalar half_angle_deg = 30.0;
        pp.query_default("fin_half_angle_deg", half_angle_deg, 30.0);
        value.half_angle = half_angle_deg * Set::Constant::Pi / 180.0;

        // -------------------------------------------------------------------
        // Station data: read as parallel arrays
        // Each station specifies: z, R1, R3, R4, is_star
        // -------------------------------------------------------------------
        std::vector<Set::Scalar> z_arr, R1_arr, R3_arr, R4_arr;
        std::vector<int>         star_arr;

        pp.queryarr("station.z",    z_arr);
        pp.queryarr("station.R1",   R1_arr);
        pp.queryarr("station.R3",   R3_arr);
        pp.queryarr("station.R4",   R4_arr);
        pp.queryarr("station.star", star_arr);

        Util::Assert(INFO, TEST(z_arr.size() == R1_arr.size()), 
                     "StarAftGrain: station.z and station.R1 must have the same length");
        Util::Assert(INFO, TEST(z_arr.size() == R3_arr.size()),
                     "StarAftGrain: station.z and station.R3 must have the same length");
        Util::Assert(INFO, TEST(z_arr.size() == R4_arr.size()),
                     "StarAftGrain: station.z and station.R4 must have the same length");
        Util::Assert(INFO, TEST(z_arr.size() == star_arr.size()),
                     "StarAftGrain: station.z and station.star must have the same length");

        value.stations.clear();
        for (int s = 0; s < (int)z_arr.size(); s++)
        {
            Station st;
            st.z       = z_arr[s];
            st.R1      = R1_arr[s];
            st.R3      = R3_arr[s];
            st.R4      = R4_arr[s];
            st.is_star = (star_arr[s] != 0);
            value.stations.push_back(st);
        }

        // Sort stations by axial position just in case
        std::sort(value.stations.begin(), value.stations.end(),
                  [](const Station& a, const Station& b){ return a.z < b.z; });

        Util::Message(INFO, "StarAftGrain: loaded ", value.stations.size(), " stations");
        Util::Message(INFO, "StarAftGrain: nfin = ", value.nfin,
                      ", half_angle = ", half_angle_deg, " deg");
        Util::Message(INFO, "StarAftGrain: R_case = ", value.R_case,
                      " m, eps = ", value.eps, " m");
    }

private:
    std::vector<Station> stations;
    Set::Scalar          eps          = 0.001;   // Interface width [m]
    Set::Scalar          R_case       = 0.127;   // Outer case radius [m]
    Set::Scalar          z_motor_lo   = 0.0;     // Head-end axial position [m]
    int                  nfin         = 6;        // Number of star fins
    Set::Scalar          half_angle   = Set::Constant::Pi / 6.0; // 30 deg in radians
};

} // namespace IC
#endif
