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
    std::array<double, NSPECIES_MAX> rho;
    int nspecies = 0;
    Set::Scalar M_normal = NAN;
    Set::Scalar M_tangent = NAN;
    Set::Scalar E = NAN;

    // Default constructor
    State() { rho.fill(0.0); }

    // Move and copy constructors
    State(const State&) = default;
    State(State&&) = default;
    State& operator=(const State&) = default;
    State& operator=(State&&) = default;

    // Constructors for convenience
    State(std::array<double, NSPECIES_MAX> a_rho, int a_nspecies, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E)
        : rho(a_rho), nspecies(a_nspecies), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E) {}

    State(Set::Patch<const Set::Scalar>& density_mf, Set::Patch<const Set::Scalar> momentum_mf, Set::Patch<const Set::Scalar> energy_mf, int i, int j, int k, int direction)
    {
        nspecies = density_mf.nComp();
        rho.fill(0.0);
        for (int n=0; n<nspecies; ++n) { rho[n] = density_mf(i,j,k,n); }
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
        else
        {
            Util::Abort(INFO, "Not supported yet");
        }
        E = energy_mf(i,j,k);
    }


    friend std::ostream &operator<<(std::ostream &os, const State &state) 
    {
        os << "nspecies=" << state.nspecies << ", ";
        for (int n=0; n<state.nspecies; ++n)
        {
            os << "rho[" << n << "]=" << state.rho[n] << ", ";
        }
        os << "Mn=" << state.M_normal << ", ";
        os << "Mt=" << state.M_tangent << ", ";
        os << "E=" << state.E << ", ";
        // do stuf
        return os;
    }
    void operator += (const State &a)
    {
        M_normal += a.M_normal; M_tangent += a.M_tangent; E += a.E;
        for (int n=0; n<nspecies; ++n) { rho[n] += a.rho[n]; }
    };
    void operator -= (const State &a)
    {
        M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E;
        for (int n=0; n<nspecies; ++n) { rho[n] -= a.rho[n]; }
    };
    void operator *= (const Set::Scalar alpha)
    {
        M_normal *= alpha; M_tangent *= alpha; E *= alpha; 
        for (int n=0; n<nspecies; ++n) { rho[n] *= alpha; }
    };
    void operator /= (const Set::Scalar alpha)
    {
        M_normal /= alpha; M_tangent /= alpha; E /= alpha; 
        for (int n=0; n<nspecies; ++n) { rho[n] /= alpha; }
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
    std::array<double, NSPECIES_MAX> mass;
    int nspecies = 0;
    Set::Scalar momentum_normal = 0.0;
    Set::Scalar momentum_tangent = 0.0;
    Set::Scalar energy = 0.0;

    // Default constructor
    Flux() { mass.fill(0.0); }

    // Move and copy constructors
    Flux(const Flux&) = default;
    Flux(Flux&&) = default;
    Flux& operator=(const Flux&) = default;
    Flux& operator=(Flux&&) = default;

    // Constructors for convenience
    Flux(std::array<double, NSPECIES_MAX> a_mass, int a_nspecies, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy) :
        mass(a_mass), nspecies(a_nspecies), momentum_normal(a_momentum_normal),
        momentum_tangent(a_momentum_tangent), energy(a_energy) {}

    Flux(const State &in) :
        mass(in.rho),
        nspecies(in.nspecies),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
        energy(in.E) {}

    Flux(State &&in) :
        mass(in.rho),
        nspecies(in.nspecies),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
        energy(in.E) {}

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
    {
        os << "nspecies=" << flux.nspecies << ", ";
        for (int n=0; n<flux.nspecies; ++n) {
            os << "mass[" << n << "]=" << flux.mass[n] << ", ";
        }
        os << "Mn=" << flux.momentum_normal << ", ";
        os << "Mt=" << flux.momentum_tangent << ", ";
        os << "E=" << flux.energy << ", ";
        // do stuff
        return os;
    }
    void operator += (const Flux &a)
    {
        for (int n=0; n<nspecies; ++n) { mass[n] += a.mass[n]; }
        momentum_normal += a.momentum_normal; momentum_tangent += a.momentum_tangent; 
        energy += a.energy;
    }
    void operator -= (const Flux &a)
    {
        for (int n=0; n<nspecies; ++n) { mass[n] -= a.mass[n]; }
        momentum_normal -= a.momentum_normal; momentum_tangent -= a.momentum_tangent; 
        energy -= a.energy;
    }
    void operator *= (const Set::Scalar alpha)
    {
        for (int n=0; n<nspecies; ++n) { mass[n] *= alpha; }
        momentum_normal *= alpha;
        momentum_tangent *= alpha;
        energy *= alpha;
    }
    void operator /= (const Set::Scalar alpha)
    {
        for (int n=0; n<nspecies; ++n) { mass[n] /= alpha; }
        momentum_normal /= alpha;
        momentum_tangent /= alpha;
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


class Riemann 
{
public:
    virtual Flux Solve(const State& lo, const State& hi, Model::Gas::Gas& gas, Set::Patch<const Set::Scalar>& X, int i, int j, int k, int side, Set::Scalar small) = 0;
    virtual ~Riemann() = default;
};


}
}
}




#endif
