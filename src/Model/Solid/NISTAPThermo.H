#ifndef MODEL_SOLID_NIST_AP_THERMO_H
#define MODEL_SOLID_NIST_AP_THERMO_H

#include <cmath>
#include "AMReX_Extension.H"
#include "AMReX_GpuQualifiers.H"

namespace Model::Solid
{
// NIST Chemistry WebBook, CAS 7790-98-9, Chase (1998).
// https://webbook.nist.gov/cgi/cbook.cgi?ID=C7790989&Type=JANAFS
// Cp: J/mol/K, Shomate H: kJ/mol, molecular weight: 117.489 g/mol.
// The two H polynomials include the crystal transition. A C1 enthalpy
// bridge preserves both endpoint enthalpies over a configurable narrow band.
// This is a diagnostic equilibrium regularization, not a transition-kinetics
// model. Constant-cp linear tails outside [298,1500] K are numerical guards.
struct NISTAPThermo
{
    static constexpr double molecular_weight = 0.117489; // kg/mol
    static constexpr double transition_temperature = 513.15;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static double PolynomialCp(double T, bool high)
    {
        const double t=T/1000.0;
        if (high)
            return (225.2150-203.0560*t+402.2310*t*t-
                    137.7330*t*t*t-7.484290/(t*t))/molecular_weight;
        return (49.81390+297.3000*t-186.1620*t*t+
                128.3380*t*t*t+0.245965/(t*t))/molecular_weight;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static double PolynomialH(double T, bool high)
    {
        const double t=T/1000.0;
        double value;
        if (high)
            value=225.2150*t-203.0560*t*t/2+402.2310*t*t*t/3-
                137.7330*t*t*t*t/4+7.484290/t-373.0730+295.7670;
        else
            value=49.81390*t+297.3000*t*t/2-186.1620*t*t*t/3+
                128.3380*t*t*t*t/4-0.245965/t-321.6170+295.7670;
        return 1000.0*value/molecular_weight;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static double H(double T, double width)
    {
        const double reference=PolynomialH(300.0,false);
        if (T<298.0)
            return PolynomialH(298.0,false)-reference+
                PolynomialCp(298.0,false)*(T-298.0);
        if (T>1500.0)
            return PolynomialH(1500.0,true)-reference+
                PolynomialCp(1500.0,true)*(T-1500.0);
        const double lo=transition_temperature-0.5*width;
        const double hi=transition_temperature+0.5*width;
        if (T<=lo) return PolynomialH(T,false)-reference;
        if (T>=hi) return PolynomialH(T,true)-reference;
        const double x=(T-lo)/width;
        const double h0=PolynomialH(lo,false), h1=PolynomialH(hi,true);
        const double c0=PolynomialCp(lo,false), c1=PolynomialCp(hi,true);
        return (2*x*x*x-3*x*x+1)*h0+(x*x*x-2*x*x+x)*width*c0+
            (-2*x*x*x+3*x*x)*h1+(x*x*x-x*x)*width*c1-reference;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static double Cp(double T, double width)
    {
        if (T<298.0) return PolynomialCp(298.0,false);
        if (T>1500.0) return PolynomialCp(1500.0,true);
        const double lo=transition_temperature-0.5*width;
        const double hi=transition_temperature+0.5*width;
        if (T<=lo) return PolynomialCp(T,false);
        if (T>=hi) return PolynomialCp(T,true);
        const double x=(T-lo)/width;
        const double h0=PolynomialH(lo,false), h1=PolynomialH(hi,true);
        const double c0=PolynomialCp(lo,false), c1=PolynomialCp(hi,true);
        return (6*x*x-6*x)*(h0-h1)/width+
            (3*x*x-4*x+1)*c0+(3*x*x-2*x)*c1;
    }
};

// Frozen-composition thermal map: all non-NIST components have common or
// individual constant cp. H is sensible + crystal-transition enthalpy,
// referenced to 300 K. Reaction Q retains the existing local-heat convention.
struct NISTAPCaloricState
{
    double solid_density=0.0;
    double other_capacity=0.0;
    double transition_width=2.0;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool Enabled() const { return solid_density>0.0; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double H(double T) const
    {
        return other_capacity*(T-300.0)+
            solid_density*NISTAPThermo::H(T,transition_width);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double Cp(double T) const
    {
        return other_capacity+
            solid_density*NISTAPThermo::Cp(T,transition_width);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double Temperature(double old_T, double energy) const
    {
        if (energy==0.0) return old_T;
        const double target=H(old_T)+energy;
        double lo=1.0, hi=5000.0;
        if (target<H(lo) || target>H(hi)) return NAN;
        double T=old_T+energy/Cp(old_T);
        if (!(T>lo && T<hi)) T=0.5*(lo+hi);
        for (int iteration=0; iteration<64; ++iteration)
        {
            const double residual=H(T)-target;
            if (std::abs(residual)<1.e-10*(1.0+std::abs(target))) return T;
            if (residual>0.0) hi=T;
            else lo=T;
            const double next=T-residual/Cp(T);
            T=(next>lo && next<hi) ? next : 0.5*(lo+hi);
        }
        return NAN; // Never silently accept an unconverged inverse.
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double DeltaT(double old_T, double energy, double fallback_capacity) const
    {
        return Enabled() ? Temperature(old_T,energy)-old_T :
            energy/fallback_capacity;
    }
};
}
#endif
