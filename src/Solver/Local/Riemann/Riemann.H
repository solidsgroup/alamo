#ifndef SOLVER_LOCAL_RIEMANN_H
#define SOLVER_LOCAL_RIEMANN_H

#include "Set/Base.H"
#include "Set/Set.H"

namespace Solver
{
namespace Local
{
namespace Riemann
{
struct State {
    Set::Scalar rho = NAN;
    Set::Scalar M_normal = NAN;
    Set::Scalar M_tangent = NAN;
    Set::Scalar E = NAN;
    std::vector<double> Y(nspecies);
    // Construtor for convenience
    State() { rho = 0.0; M_normal = 0.0, M_tangent = 0.0; E = 0.0; Y = std::vector<double>(nspecies, 0.0);}
    State(Set::Scalar a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E, std::vector<double> a_Y)
        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E), Y(a_Y) {}
    State(Set::Patch<const Set::Scalar> density_mf, Set::Patch<const Set::Scalar> momentum_mf, Set::Patch<const Set::Scalar> energy_mf, Set::Patch<const Set::Scalar> massfraction_mf, 
          int i, int j, int k, int direction)
    {
        rho = density_mf(i,j,k);
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
        for (int ii=0; ii<nspecies; ++ii) {
            Y[ii] = massfraction_mf(i,j,k,ii);
        }
    }

    friend std::ostream &operator<<(std::ostream &os, const State &state) 
    {
        os << "rho=" << state.rho << ", "; 
        os << "Mn=" << state.M_normal << ", ";
        os << "Mt=" << state.M_tangent << ", ";
        os << "E=" << state.E << ", ";
        for (int ii=0; ii<nspecies; ++ii) {
            os << "Y[" << ii << "]=" << state.Y[ii] << ", ";
        }
        // do stuf
        return os;
    }
    void operator += (const State &a)
    {
        rho += a.rho; M_normal += a.M_normal; M_tangent += a.M_tangent; E += a.E;
        for (int ii=0; ii<nspecies; ++ii) { Y[ii] += a.Y[ii]; }
    };
    void operator -= (const State &a)
    {
        rho -= a.rho; M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E; 
        for (int ii=0; ii<nspecies; ++ii) { Y[ii] -= a.Y[ii]; }
    };
    void operator *= (const Set::Scalar alpha)
    {
        rho *= alpha; M_normal *= alpha; M_tangent *= alpha; E *= alpha; 
        for (int ii=0; ii<nspecies; ++ii) { Y[ii] *= alpha; }
    };
    void operator /= (const Set::Scalar alpha)
    {
        rho /= alpha; M_normal /= alpha; M_tangent /= alpha; E /= alpha; 
        for (int ii=0; ii<nspecies; ++ii) { Y[ii] /= alpha; }
    };
    friend State operator + (const State &a, const State &b)
    {
        return State(a.rho+b.rho, a.M_normal+b.M_normal, a.M_tangent+b.M_tangent, a.E+b.E, std::transform(a.Y.begin(), a.Y.end(), b.Y.begin(), b.Y.begin(), std::plus<double>()));
    }
    friend State operator - (const State &a, const State &b)
    {
        return State(a.rho-b.rho, a.M_normal-b.M_normal, a.M_tangent-b.M_tangent, a.E-b.E, std::transform(a.Y.begin(), a.Y.end(), b.Y.begin(), b.Y.begin(), std::minus<double>()));
    }
    friend State operator * (const Set::Scalar alpha, const State &b)
    {
        return State(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha, std::transform(b.Y.begin(), b.Y.end(), b.Y.begin(), [alpha](double val){ return val*alpha; }));
    }
    friend State operator * (const State &b, const Set::Scalar alpha)
    {
        return State(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha, std::transform(b.Y.begin(), b.Y.end(), b.Y.begin(), [alpha](double val){ return val*alpha; }));
    }
    friend State operator / (const State &b, const Set::Scalar alpha)
    {
        return State(b.rho/alpha, b.M_normal/alpha, b.M_tangent/alpha, b.E/alpha, std::transform(b.Y.begin(), b.Y.end(), b.Y.begin(), [alpha](double val){ return val/alpha; }));
    }
};

struct Flux {
    Set::Scalar mass;
    Set::Scalar momentum_normal;
    Set::Scalar momentum_tangent;
    Set::Scalar energy;
    std::vector<double> rhoY(nspecies);
    Flux() : mass(0.0), momentum_normal(0.0), momentum_tangent(0.0), energy(0.0),  {}
    Flux(Set::Scalar a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy) :
        mass(a_mass), momentum_normal(a_momentum_normal),
        momentum_tangent(a_momentum_tangent), energy(a_energy) {}

    Flux(const State &in) :
        mass(in.rho),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
        energy(in.E) {}

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
    {
        os << "mass=" << flux.mass << ", "; 
        os << "Mn=" << flux.momentum_normal << ", ";
        os << "Mt=" << flux.momentum_tangent << ", ";
        os << "E=" << flux.energy << ", ";
        // do stuff
        return os;
    }
    void operator += (const Flux &a)
    {
        mass += a.mass;
        momentum_normal += a.momentum_normal; momentum_tangent += a.momentum_tangent; 
        energy += a.energy;
    }
    void operator -= (const Flux &a)
    {
        mass -= a.mass;
        momentum_normal -= a.momentum_normal; momentum_tangent -= a.momentum_tangent; 
        energy -= a.energy;
    }
    void operator *= (const Set::Scalar alpha)
    {
        mass *= alpha;
        momentum_normal *= alpha;
        momentum_tangent *= alpha;
        energy *= alpha;
    }
    void operator /= (const Set::Scalar alpha)
    {
        mass /= alpha;
        momentum_normal /= alpha;
        momentum_tangent /= alpha;
        energy /= alpha;
    }
    friend Flux operator + (const Flux &a, const Flux &b)
    {
        return Flux(a.mass+b.mass, a.momentum_normal+b.momentum_normal,
                    a.momentum_tangent+b.momentum_tangent, a.energy+b.energy);
    }
    friend Flux operator - (const Flux &a, const Flux &b)
    {
        return Flux(a.mass-b.mass, a.momentum_normal-b.momentum_normal,
                    a.momentum_tangent-b.momentum_tangent, a.energy-b.energy);
    }
    friend Flux operator * (const Flux &a, const Set::Scalar beta)
    {
        return Flux(a.mass*beta, a.momentum_normal*beta,
                    a.momentum_tangent*beta, a.energy*beta);
    }
    friend Flux operator * (const Set::Scalar beta,const Flux &a)
    {
        return Flux(a.mass*beta, a.momentum_normal*beta,
                    a.momentum_tangent*beta, a.energy*beta);
    }
    friend Flux operator / (const Flux &a, const Set::Scalar beta)
    {
        return Flux(a.mass/beta, a.momentum_normal/beta,
                    a.momentum_tangent/beta, a.energy/beta);
    }

};


class Riemann 
{
public:
    virtual Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small) = 0;
    virtual ~Riemann() = default;
};


}
}
}




#endif
