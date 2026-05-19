#ifndef MODEL_CHEMISTRY_ROCFIRE_H
#define MODEL_CHEMISTRY_ROCFIRE_H

#include "AMReX_Extension.H"
#include "Model/Chemistry/Chemistry.H"
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace Model
{
namespace Chemistry
{

struct Rocfire
{
    static constexpr const char* name = "rocfire";
    static constexpr int nreactions = 4;
    static constexpr int AP = 0;
    static constexpr int Binder = 1;
    static constexpr int Mono = 2;
    static constexpr int Premixed = 3;
    static constexpr int Primary = 4;
    static constexpr int Final = 5;

    std::array<double, nreactions> A = {{27000.0, 25300.0, 5000.0, 400.0}};
    std::array<double, nreactions> pressure_exponent = {{1.6, 1.7, 1.7, 1.7}};
    std::array<double, nreactions> activation_energy =
        {{15.9 * 4184.0, 14.5 * 4184.0, 11.3 * 4184.0, 14.5 * 4184.0}};
    std::array<double, nreactions> qgas =
        {{430.0 * 4184.0, 447.0 * 4184.0, 8365.0 * 4184.0, 2127.0 * 4184.0}};
    std::array<double, 2> qsolid = {{-100.0 * 4184.0, -300.0 * 4184.0}};
    double pressure_reference = 100000.0;
    int nsubsteps = 1;

    Rocfire() = default;
    Rocfire(IO::ParmParse &pp, std::string name)
    {
        pp_queryclass(name,*this);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    static double sum_species(const std::array<double, NSPECIES>& values)
    {
        double sum = 0.0;
        for (int n = 0; n < NSPECIES; ++n) sum += values[n];
        return sum;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::array<double, nreactions>
    compute_progress_rates(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY) const
    {
        std::array<double, nreactions> progress = {{0.0, 0.0, 0.0, 0.0}};

        double rho = 0.0;
        std::array<double, NSPECIES> Y;
        for (int n = 0; n < NSPECIES; ++n)
        {
            rho += rhoY[n];
        }
        if (!(rho > 0.0) || !(P > 0.0) || !(T > 0.0) || !(pressure_reference > 0.0))
        {
            return progress;
        }

        for (int n = 0; n < NSPECIES; ++n)
        {
            Y[n] = std::max(rhoY[n], 0.0) / rho;
        }

        std::array<double, nreactions> k;
        for (int r = 0; r < nreactions; ++r)
        {
            k[r] = A[r] * std::pow(P / pressure_reference, pressure_exponent[r])
                 * std::exp(-activation_energy[r] / (Set::Constant::Rg * T));
        }

        progress[0] = k[0] * Y[AP];
        progress[1] = k[1] * Y[Binder];
        progress[2] = k[2] * Y[AP] * Y[Binder];
        progress[3] = k[3] * Y[Mono] * Y[Premixed];
        return progress;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double, NSPECIES>, Set::Scalar>
    source_from_progress(const std::array<double, nreactions>& progress) const
    {
        std::array<double, NSPECIES> wdot;
        for (int n = 0; n < NSPECIES; ++n) wdot[n] = 0.0;

        wdot[AP]       -= progress[0];
        wdot[Mono]     += progress[0];

        wdot[Binder]   -= progress[1];
        wdot[Premixed] += progress[1];

        wdot[AP]       -= 7.974 * progress[2];
        wdot[Binder]   -= progress[2];
        wdot[Primary]  += 8.974 * progress[2];

        wdot[Mono]     -= 3.19 * progress[3];
        wdot[Premixed] -= progress[3];
        wdot[Final]    += 4.19 * progress[3];

        const double qdot =
            (qgas[0] + qsolid[0]) * progress[0]
          + (qgas[1] + qsolid[1]) * progress[1]
          + (8.974 * qgas[2] + 7.974 * qsolid[0] + qsolid[1]) * progress[2]
          + (4.19 * qgas[3]) * progress[3];

        return {wdot, qdot};
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double, NSPECIES>, Set::Scalar>
    evaluate_source(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt) const
    {
        std::array<double, nreactions> progress = compute_progress_rates(P, T, rhoY);
        std::pair<std::array<double, NSPECIES>, Set::Scalar> source = source_from_progress(progress);

        if (dt > 0.0)
        {
            double limiter = 1.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                if (source.first[n] < 0.0)
                {
                    const double available = std::max(rhoY[n], 0.0);
                    const double requested = -source.first[n] * dt;
                    if (requested > available && requested > 0.0)
                    {
                        limiter = std::min(limiter, available / requested);
                    }
                }
            }

            if (limiter < 1.0)
            {
                for (int n = 0; n < NSPECIES; ++n) source.first[n] *= limiter;
                source.second *= limiter;
            }
        }

        return source;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    std::pair<std::array<double,NSPECIES>,Set::Scalar>
    compute(
        const Set::Scalar P,
        const Set::Scalar T,
        const std::array<double, NSPECIES>& rhoY,
        const Set::Scalar dt,
        const Model::Gas::Gas* /*gas*/) const
    {
        std::array<double, NSPECIES> wdot;
        for (int n = 0; n < NSPECIES; ++n)
        {
            wdot[n] = 0.0;
        }

        if constexpr (NSPECIES != 6)
        {
            Util::Abort(INFO, "Rocfire chemistry requires NSPECIES=6.");
            return {wdot, 0.0};
        }
        else
        {
            std::array<double, NSPECIES> rhoY_valid;
            double density_scale = 0.0;
            for (int n = 0; n < NSPECIES; ++n)
            {
                if (!std::isfinite(rhoY[n]))
                {
                    Util::Abort(INFO, "Rocfire chemistry received non-finite species density for species ",
                        n, ": ", rhoY[n]);
                }
                density_scale += std::abs(rhoY[n]);
            }

            const double density_tol = 1.0e-12 * ((density_scale > 1.0) ? density_scale : 1.0);
            for (int n = 0; n < NSPECIES; ++n)
            {
                rhoY_valid[n] = rhoY[n];
                if (rhoY_valid[n] < 0.0)
                {
                    if (rhoY_valid[n] < -density_tol)
                    {
                        Util::Abort(INFO, "Rocfire chemistry received negative species density for species ",
                            n, ": rhoY=", rhoY_valid[n], ", tolerance=", density_tol);
                    }
                    rhoY_valid[n] = 0.0;
                }
            }

            if (!(sum_species(rhoY_valid) > 0.0)) return {wdot, 0.0};

            if (dt > 0.0 && nsubsteps > 1)
            {
                const double dt_sub = dt / nsubsteps;
                std::array<double, NSPECIES> rhoY_sub = rhoY_valid;
                Set::Scalar qdot_integral = 0.0;

                for (int s = 0; s < nsubsteps; ++s)
                {
                    std::pair<std::array<double, NSPECIES>, Set::Scalar> source =
                        evaluate_source(P, T, rhoY_sub, dt_sub);
                    for (int n = 0; n < NSPECIES; ++n)
                    {
                        rhoY_sub[n] += source.first[n] * dt_sub;
                        if (rhoY_sub[n] < 0.0) rhoY_sub[n] = 0.0;
                    }
                    qdot_integral += source.second * dt_sub;
                }

                for (int n = 0; n < NSPECIES; ++n)
                {
                    wdot[n] = (rhoY_sub[n] - rhoY_valid[n]) / dt;
                }

                return {wdot, qdot_integral / dt};
            }

            return evaluate_source(P, T, rhoY_valid, dt);
        }
    }

    static void Parse(Rocfire& value, IO::ParmParse& pp)
    {
        Util::Message(INFO, pp.getPrefix());
        Util::Assert(INFO, TEST(NSPECIES == 6), "Rocfire chemistry requires NSPECIES=6.");

        std::vector<double> A;
        std::vector<double> pressure_exponent;
        std::vector<double> activation_energy;
        std::vector<double> qgas;
        std::vector<double> qsolid;

        pp.queryarr_default("A", A,
            {"27.0_g/cm^3/s", "25.3_g/cm^3/s", "5.0_g/cm^3/s", "0.4_g/cm^3/s"},
            Unit::Density() / Unit::Time());
        pp.queryarr_default("pressure_exponent", pressure_exponent, {1.6, 1.7, 1.7, 1.7});
        pp.queryarr_default("activation_energy", activation_energy,
            {"15.9_kcal/mol", "14.5_kcal/mol", "11.3_kcal/mol", "14.5_kcal/mol"},
            Unit::Energy() / Unit::Amount());
        pp.queryarr_default("qgas", qgas,
            {"430.0_cal/g", "447.0_cal/g", "8365.0_cal/g", "2127.0_cal/g"},
            Unit::Energy() / Unit::Mass());
        pp.queryarr_default("qsolid", qsolid,
            {"-100.0_cal/g", "-300.0_cal/g"},
            Unit::Energy() / Unit::Mass());
        pp.query_default("pressure_reference", value.pressure_reference, "1_bar", Unit::Pressure());
        pp.query_default("nsubsteps", value.nsubsteps, 1);

        Util::Assert(INFO, TEST(A.size() == nreactions), "chemistry.rocfire.A must have 4 entries.");
        Util::Assert(INFO, TEST(pressure_exponent.size() == nreactions),
            "chemistry.rocfire.pressure_exponent must have 4 entries.");
        Util::Assert(INFO, TEST(activation_energy.size() == nreactions),
            "chemistry.rocfire.activation_energy must have 4 entries.");
        Util::Assert(INFO, TEST(qgas.size() == nreactions), "chemistry.rocfire.qgas must have 4 entries.");
        Util::Assert(INFO, TEST(qsolid.size() == 2), "chemistry.rocfire.qsolid must have 2 entries.");
        Util::Assert(INFO, TEST(value.pressure_reference > 0.0),
            "chemistry.rocfire.pressure_reference must be positive.");
        Util::Assert(INFO, TEST(value.nsubsteps > 0), "chemistry.rocfire.nsubsteps must be > 0.");

        std::copy(A.begin(), A.end(), value.A.begin());
        std::copy(pressure_exponent.begin(), pressure_exponent.end(), value.pressure_exponent.begin());
        std::copy(activation_energy.begin(), activation_energy.end(), value.activation_energy.begin());
        std::copy(qgas.begin(), qgas.end(), value.qgas.begin());
        std::copy(qsolid.begin(), qsolid.end(), value.qsolid.begin());
    }

};

} // namespace Chemistry
} // namespace Model

#endif
