#ifndef SOLVER_LOCAL_RIEMANN_H
#define SOLVER_LOCAL_RIEMANN_H

#include "Set/Base.H"
#include "Set/Set.H"
#include "Model/Gas/Gas.H"
#include <memory>
#include <array>

namespace Solver
{
namespace Local
{
namespace Riemann
{
struct State {
    std::array<double, NSPECIES> rho;
    Set::Scalar M_normal = NAN;
    Set::Scalar M_tangent = NAN;
#if AMREX_SPACEDIM == 3
    Set::Scalar M_tangent2 = NAN;
#endif
    Set::Scalar E = NAN;

    // Default constructor
    State() { rho.fill(0.0); }

    // Move and copy constructors
    State(const State&) = default;
    State(State&&) = default;
    State& operator=(const State&) = default;
    State& operator=(State&&) = default;

    // Constructors for convenience
    State(std::array<double, NSPECIES> a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E)
        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E) {}

#if AMREX_SPACEDIM == 3
    State(std::array<double, NSPECIES> a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_M_tangent2, Set::Scalar a_E)
        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), M_tangent2(a_M_tangent2), E(a_E) {}
#endif

    State(Set::Patch<const Set::Scalar>& density_mf, Set::Patch<const Set::Scalar> momentum_mf, Set::Patch<const Set::Scalar> energy_mf, int i, int j, int k, int direction)
    {
        rho.fill(0.0);
        for (int n=0; n<NSPECIES; ++n) { rho[n] = density_mf(i,j,k,n); }
        if (direction == 0)
        {
            M_normal = momentum_mf(i,j,k,0);
            M_tangent = momentum_mf(i,j,k,1);
        }
        else if (direction == 1)
        {
            M_normal = momentum_mf(i,j,k,1);
            M_tangent = momentum_mf(i,j,k,0);
        }
#if AMREX_SPACEDIM == 3
        else if (direction == 2)
        {
            M_normal = momentum_mf(i,j,k,2);
            M_tangent = momentum_mf(i,j,k,0);
            M_tangent2 = momentum_mf(i,j,k,1);
        }
#endif
        else
        {
            Util::Abort(INFO, "Not supported yet");
        }
#if AMREX_SPACEDIM == 3
        if (direction == 0 || direction == 1) M_tangent2 = momentum_mf(i,j,k,2);
#endif
        E = energy_mf(i,j,k);
    }


    friend std::ostream &operator<<(std::ostream &os, const State &state) 
    {
        os << "NSPECIES=" << NSPECIES << ", ";
        for (int n=0; n<NSPECIES; ++n)
        {
            os << "rho[" << n << "]=" << state.rho[n] << ", ";
        }
        os << "Mn=" << state.M_normal << ", ";
        os << "Mt=" << state.M_tangent << ", ";
#if AMREX_SPACEDIM == 3
        os << "Mt2=" << state.M_tangent2 << ", ";
#endif
        os << "E=" << state.E << ", ";
        // do stuf
        return os;
    }
    void operator += (const State &a)
    {
        M_normal += a.M_normal; M_tangent += a.M_tangent; E += a.E;
#if AMREX_SPACEDIM == 3
        M_tangent2 += a.M_tangent2;
#endif
        for (int n=0; n<NSPECIES; ++n) { rho[n] += a.rho[n]; }
    };
    void operator -= (const State &a)
    {
        M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E;
#if AMREX_SPACEDIM == 3
        M_tangent2 -= a.M_tangent2;
#endif
        for (int n=0; n<NSPECIES; ++n) { rho[n] -= a.rho[n]; }
    };
    void operator *= (const Set::Scalar alpha)
    {
        M_normal *= alpha; M_tangent *= alpha; E *= alpha; 
#if AMREX_SPACEDIM == 3
        M_tangent2 *= alpha;
#endif
        for (int n=0; n<NSPECIES; ++n) { rho[n] *= alpha; }
    };
    void operator /= (const Set::Scalar alpha)
    {
        M_normal /= alpha; M_tangent /= alpha; E /= alpha; 
#if AMREX_SPACEDIM == 3
        M_tangent2 /= alpha;
#endif
        for (int n=0; n<NSPECIES; ++n) { rho[n] /= alpha; }
    };
    friend State operator + (State a, const State &b)
    {
        a += b;
        return a;
    }
    friend State operator - (State a, const State &b)
    {
        a -= b;
        return a;
    }
    friend State operator * (const Set::Scalar alpha, State b)
    {
        b *= alpha;
        return b;
    }
    friend State operator * (State b, const Set::Scalar alpha)
    {
        b *= alpha;
        return b;
    }
    friend State operator / (State b, const Set::Scalar alpha)
    {
        b /= alpha;
        return b;
    }
};

struct Flux {
    std::array<double, NSPECIES> mass;
    Set::Scalar momentum_normal = 0.0;
    Set::Scalar momentum_tangent = 0.0;
#if AMREX_SPACEDIM == 3
    Set::Scalar momentum_tangent2 = 0.0;
#endif
    Set::Scalar energy = 0.0;

    // Default constructor
    Flux() { mass.fill(0.0); }

    // Move and copy constructors
    Flux(const Flux&) = default;
    Flux(Flux&&) = default;
    Flux& operator=(const Flux&) = default;
    Flux& operator=(Flux&&) = default;

    // Constructors for convenience
    Flux(std::array<double, NSPECIES> a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy) :
        mass(a_mass), momentum_normal(a_momentum_normal),
        momentum_tangent(a_momentum_tangent), energy(a_energy) {}

#if AMREX_SPACEDIM == 3
    Flux(std::array<double, NSPECIES> a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_momentum_tangent2, Set::Scalar a_energy) :
        mass(a_mass), momentum_normal(a_momentum_normal),
        momentum_tangent(a_momentum_tangent), momentum_tangent2(a_momentum_tangent2), energy(a_energy) {}
#endif

    Flux(const State &in) :
        mass(in.rho),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
#if AMREX_SPACEDIM == 3
        momentum_tangent2(in.M_tangent2),
#endif
        energy(in.E) {}

    Flux(State &&in) :
        mass(in.rho),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
#if AMREX_SPACEDIM == 3
        momentum_tangent2(in.M_tangent2),
#endif
        energy(in.E) {}

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
    {
        os << "NSPECIES=" << NSPECIES << ", ";
        for (int n=0; n<NSPECIES; ++n) {
            os << "mass[" << n << "]=" << flux.mass[n] << ", ";
        }
        os << "Mn=" << flux.momentum_normal << ", ";
        os << "Mt=" << flux.momentum_tangent << ", ";
#if AMREX_SPACEDIM == 3
        os << "Mt2=" << flux.momentum_tangent2 << ", ";
#endif
        os << "E=" << flux.energy << ", ";
        // do stuff
        return os;
    }
    void operator += (const Flux &a)
    {
        for (int n=0; n<NSPECIES; ++n) { mass[n] += a.mass[n]; }
        momentum_normal += a.momentum_normal; momentum_tangent += a.momentum_tangent; 
#if AMREX_SPACEDIM == 3
        momentum_tangent2 += a.momentum_tangent2;
#endif
        energy += a.energy;
    }
    void operator -= (const Flux &a)
    {
        for (int n=0; n<NSPECIES; ++n) { mass[n] -= a.mass[n]; }
        momentum_normal -= a.momentum_normal; momentum_tangent -= a.momentum_tangent; 
#if AMREX_SPACEDIM == 3
        momentum_tangent2 -= a.momentum_tangent2;
#endif
        energy -= a.energy;
    }
    void operator *= (const Set::Scalar alpha)
    {
        for (int n=0; n<NSPECIES; ++n) { mass[n] *= alpha; }
        momentum_normal *= alpha;
        momentum_tangent *= alpha;
#if AMREX_SPACEDIM == 3
        momentum_tangent2 *= alpha;
#endif
        energy *= alpha;
    }
    void operator /= (const Set::Scalar alpha)
    {
        for (int n=0; n<NSPECIES; ++n) { mass[n] /= alpha; }
        momentum_normal /= alpha;
        momentum_tangent /= alpha;
#if AMREX_SPACEDIM == 3
        momentum_tangent2 /= alpha;
#endif
        energy /= alpha;
    }
    friend Flux operator + (Flux a, const Flux &b)
    {
        a += b;
        return a;
    }
    friend Flux operator - (Flux a, const Flux &b)
    {
        a -= b;
        return a;
    }
    friend Flux operator * (Flux a, const Set::Scalar beta)
    {
        a *= beta;
        return a;
    }
    friend Flux operator * (const Set::Scalar beta, Flux a)
    {
        a *= beta;
        return a;
    }
    friend Flux operator / (Flux a, const Set::Scalar beta)
    {
        a /= beta;
        return a;
    }

};

struct FluxDivergence {
    std::array<Set::Scalar, NSPECIES> mass;
    Set::Vector momentum;
    Set::Scalar energy = 0.0;

    FluxDivergence()
        : momentum(Set::Vector::Zero())
    {
        mass.fill(0.0);
    }

    FluxDivergence(
        std::array<Set::Scalar, NSPECIES> a_mass,
        Set::Vector a_momentum,
        Set::Scalar a_energy)
        : mass(a_mass), momentum(a_momentum), energy(a_energy) {}
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar ComputeStateTemperature(
    Model::Gas::Gas& gas,
    Set::Scalar density,
    const State& state,
    Set::Scalar energy,
    Set::Scalar temperature_guess,
    Set::Patch<const Set::Scalar>& X,
    int i,
    int j,
    int k)
{
#if AMREX_SPACEDIM == 2
    return gas.ComputeT(density, state.M_normal, state.M_tangent, energy, temperature_guess, X, i, j, k);
#elif AMREX_SPACEDIM == 3
    return gas.ComputeT(density, state.M_normal, state.M_tangent, state.M_tangent2, energy, temperature_guess, X, i, j, k);
#endif
}


class Riemann 
{
public:
    virtual Flux Solve(const State& lo, const State& hi, Model::Gas::Gas& gas, Set::Patch<const Set::Scalar>& X, int i, int j, int k, int side, Set::Scalar small) = 0;

    FluxDivergence ComputeFluxDivergence(
        const State& xlo, const State& x, const State& xhi,
        const State& ylo, const State& y, const State& yhi,
#if AMREX_SPACEDIM == 3
        const State& zlo, const State& z, const State& zhi,
#endif
        Model::Gas::Gas& gas,
        Set::Patch<const Set::Scalar>& X,
        int i,
        int j,
        int k,
        Set::Scalar eta,
        const Set::Scalar* DX,
        Set::Scalar small)
    {
        Flux flux_xlo = Solve(xlo, x, gas, X, i, j, k, 0, small) * eta;
        Flux flux_xhi = Solve(x, xhi, gas, X, i, j, k, 1, small) * eta;
        Flux flux_ylo = Solve(ylo, y, gas, X, i, j, k, 2, small) * eta;
        Flux flux_yhi = Solve(y, yhi, gas, X, i, j, k, 3, small) * eta;
#if AMREX_SPACEDIM == 3
        Flux flux_zlo = Solve(zlo, z, gas, X, i, j, k, 4, small) * eta;
        Flux flux_zhi = Solve(z, zhi, gas, X, i, j, k, 5, small) * eta;
#endif

        FluxDivergence result;
        for (int n=0; n<NSPECIES; ++n)
        {
            result.mass[n] =
                (flux_xlo.mass[n] - flux_xhi.mass[n]) / DX[0] +
                (flux_ylo.mass[n] - flux_yhi.mass[n]) / DX[1]
#if AMREX_SPACEDIM == 3
                + (flux_zlo.mass[n] - flux_zhi.mass[n]) / DX[2]
#endif
                ;
        }

        result.momentum(0) =
            (flux_xlo.momentum_normal  - flux_xhi.momentum_normal ) / DX[0] +
            (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1]
#if AMREX_SPACEDIM == 3
            + (flux_zlo.momentum_tangent - flux_zhi.momentum_tangent) / DX[2]
#endif
            ;

        result.momentum(1) =
            (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
            (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1]
#if AMREX_SPACEDIM == 3
            + (flux_zlo.momentum_tangent2 - flux_zhi.momentum_tangent2) / DX[2]
#endif
            ;

#if AMREX_SPACEDIM == 3
        result.momentum(2) =
            (flux_xlo.momentum_tangent2 - flux_xhi.momentum_tangent2) / DX[0] +
            (flux_ylo.momentum_tangent2 - flux_yhi.momentum_tangent2) / DX[1] +
            (flux_zlo.momentum_normal   - flux_zhi.momentum_normal  ) / DX[2];
#endif

        result.energy =
            (flux_xlo.energy - flux_xhi.energy) / DX[0] +
            (flux_ylo.energy - flux_yhi.energy) / DX[1]
#if AMREX_SPACEDIM == 3
            + (flux_zlo.energy - flux_zhi.energy) / DX[2]
#endif
            ;

        return result;
    }

    virtual ~Riemann() = default;
};


}
}
}




#endif
