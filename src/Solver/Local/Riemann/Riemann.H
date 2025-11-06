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
//struct State {
//    Set::Scalar rho = NAN;
//    Set::Scalar M_normal = NAN;
//    Set::Scalar M_tangent = NAN;
//    Set::Scalar E = NAN;
//    // Construtor for convenience
//    State() { rho = 0.0; M_normal = 0.0, M_tangent = 0.0; E = 0.0; }
//    State(Set::Scalar a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E)
//        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E) {}
//    State(Set::Patch<const Set::Scalar> density_mf,
//          Set::Patch<const Set::Scalar> momentum_mf,
//          Set::Patch<const Set::Scalar> energy_mf,
//          int i, int j, int k, int direction)
//    {
//        rho = density_mf(i,j,k);
//        if (direction == 0)
//        {
//            M_normal = momentum_mf(i,j,k,0);
//            M_tangent = momentum_mf(i,j,k,1);
//        }
//        else if (direction == 1)
//        {
//            M_normal = momentum_mf(i,j,k,1);
//            M_tangent = momentum_mf(i,j,k,0);
//        }
//        else
//        {
//            Util::Abort(INFO, "Not supported yet");
//        }
//        E = energy_mf(i,j,k);
//    }
//
//    friend std::ostream &operator<<(std::ostream &os, const State &state) 
//    {
//        os << "rho=" << state.rho << ", "; 
//        os << "Mn=" << state.M_normal << ", ";
//        os << "Mt=" << state.M_tangent << ", ";
//        os << "E=" << state.E << ", ";
//        // do stuf
//        return os;
//    }
//    void operator += (const State &a)
//    {
//        rho += a.rho; M_normal += a.M_normal; M_tangent += a.M_tangent; E += a.E;
//    };
//    void operator -= (const State &a)
//    {
//        rho -= a.rho; M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E; 
//    };
//    void operator *= (const Set::Scalar alpha)
//    {
//        rho *= alpha; M_normal *= alpha; M_tangent *= alpha; E *= alpha; 
//    };
//    void operator /= (const Set::Scalar alpha)
//    {
//        rho /= alpha; M_normal /= alpha; M_tangent /= alpha; E /= alpha; 
//    };
//    friend State operator + (const State &a, const State &b)
//    {
//        return State(a.rho+b.rho, a.M_normal+b.M_normal, a.M_tangent+b.M_tangent, a.E+b.E);
//    }
//    friend State operator - (const State &a, const State &b)
//    {
//        return State(a.rho-b.rho, a.M_normal-b.M_normal, a.M_tangent-b.M_tangent, a.E-b.E);
//    }
//    friend State operator * (const Set::Scalar alpha, const State &b)
//    {
//        return State(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha);
//    }
//    friend State operator * (const State &b, const Set::Scalar alpha)
//    {
//        return State(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha);
//    }
//    friend State operator / (const State &b, const Set::Scalar alpha)
//    {
//        return State(b.rho/alpha, b.M_normal/alpha, b.M_tangent/alpha, b.E/alpha);
//    }
//};

struct State {
    std::vector<double> rho;
    Set::Scalar M_normal = NAN;
    Set::Scalar M_tangent = NAN;
    Set::Scalar E = NAN;
    Set::Scalar T = NAN;
    // Construtor for convenience
    State() { rho = std::vector<double>(0); M_normal = 0.0, M_tangent = 0.0; E = 0.0; T = 0.0;}
    State(std::vector<double> a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E )
        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E) {}
    State(Set::Patch<const Set::Scalar>& density_mf,
          Set::Patch<const Set::Scalar> momentum_mf,
          Set::Patch<const Set::Scalar> energy_mf,
          int i, int j, int k, int direction)
    {
        size_t nspecies = density_mf.nComp();
        for (size_t n=0; n<nspecies; ++n) { rho.push_back(density_mf(i,j,k,n)); }
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
        os << "nspecies=" << state.rho.size() << ", ";
        for (size_t n=0; n<state.rho.size(); ++n) {
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
        for (size_t n=0; n<rho.size(); ++n) { rho[n] += a.rho[n]; }
    };
    void operator -= (const State &a)
    {
        M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E;
        for (size_t n=0; n<rho.size(); ++n) { rho[n] -= a.rho[n]; }
    };
    void operator *= (const Set::Scalar alpha)
    {
        M_normal *= alpha; M_tangent *= alpha; E *= alpha;
        for (size_t n=0; n<rho.size(); ++n) { rho[n] *= alpha; }
    };
    void operator /= (const Set::Scalar alpha)
    {
        M_normal /= alpha; M_tangent /= alpha; E /= alpha;
        for (size_t n=0; n<rho.size(); ++n) { rho[n] /= alpha; }
    };
    friend State operator + (const State &a, const State &b)
    {
        std::vector<double> rho_sum(a.rho.size());
        std::transform(a.rho.begin(), a.rho.end(), b.rho.begin(), rho_sum.begin(), std::plus<double>());
        return State(rho_sum, a.M_normal+b.M_normal, a.M_tangent+b.M_tangent, a.E+b.E);
    }
    friend State operator - (const State &a, const State &b)
    {
        std::vector<double> rho_diff(a.rho.size());
        std::transform(a.rho.begin(), a.rho.end(), b.rho.begin(), rho_diff.begin(), std::minus<double>());
        return State(rho_diff, a.M_normal-b.M_normal, a.M_tangent-b.M_tangent, a.E-b.E);
    }
    friend State operator * (const Set::Scalar alpha, const State &b)
    {
        std::vector<double> rho_prod(b.rho.size());
        std::transform(b.rho.begin(), b.rho.end(), rho_prod.begin(), [alpha](double val){ return val*alpha; });
        return State(rho_prod, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha);
    }
    friend State operator * (const State &b, const Set::Scalar alpha)
    {
        std::vector<double> rho_prod(b.rho.size());
        std::transform(b.rho.begin(), b.rho.end(), rho_prod.begin(), [alpha](double val){ return val*alpha; });
        return State(rho_prod, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha);
    }
    friend State operator / (const State &b, const Set::Scalar alpha)
    {
        std::vector<double> rho_quot(b.rho.size());
        std::transform(b.rho.begin(), b.rho.end(), rho_quot.begin(), [alpha](double val){ return val/alpha; });
        return State(rho_quot, b.M_normal/alpha, b.M_tangent/alpha, b.E/alpha);
    }
};

//struct Flux {
//    Set::Scalar mass;
//    Set::Scalar momentum_normal;
//    Set::Scalar momentum_tangent;
//    Set::Scalar energy;
//    Flux() : mass(0.0), momentum_normal(0.0), momentum_tangent(0.0), energy(0.0)  {}
//    Flux(Set::Scalar a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy) :
//        mass(a_mass), momentum_normal(a_momentum_normal), momentum_tangent(a_momentum_tangent), energy(a_energy) {}
//
//    Flux(const State &in) :
//        mass(in.rho),
//        momentum_normal(in.M_normal),
//        momentum_tangent(in.M_tangent),
//        energy(in.E) {}
//
//    friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
//    {
//        os << "mass=" << flux.mass << ", "; 
//        os << "Mn=" << flux.momentum_normal << ", ";
//        os << "Mt=" << flux.momentum_tangent << ", ";
//        os << "E=" << flux.energy << ", ";
//        // do stuff
//        return os;
//    }
//    void operator += (const Flux &a)
//    {
//        mass += a.mass;
//        momentum_normal += a.momentum_normal; momentum_tangent += a.momentum_tangent; 
//        energy += a.energy;
//    }
//    void operator -= (const Flux &a)
//    {
//        mass -= a.mass;
//        momentum_normal -= a.momentum_normal; momentum_tangent -= a.momentum_tangent; 
//        energy -= a.energy;
//    }
//    void operator *= (const Set::Scalar alpha)
//    {
//        mass *= alpha;
//        momentum_normal *= alpha;
//        momentum_tangent *= alpha;
//        energy *= alpha;
//    }
//    void operator /= (const Set::Scalar alpha)
//    {
//        mass /= alpha;
//        momentum_normal /= alpha;
//        momentum_tangent /= alpha;
//        energy /= alpha;
//    }
//    friend Flux operator + (const Flux &a, const Flux &b)
//    {
//        return Flux(a.mass+b.mass, a.momentum_normal+b.momentum_normal,
//                    a.momentum_tangent+b.momentum_tangent, a.energy+b.energy);
//    }
//    friend Flux operator - (const Flux &a, const Flux &b)
//    {
//        return Flux(a.mass-b.mass, a.momentum_normal-b.momentum_normal,
//                    a.momentum_tangent-b.momentum_tangent, a.energy-b.energy);
//    }
//    friend Flux operator * (const Flux &a, const Set::Scalar beta)
//    {
//        return Flux(a.mass*beta, a.momentum_normal*beta,
//                    a.momentum_tangent*beta, a.energy*beta);
//    }
//    friend Flux operator * (const Set::Scalar beta,const Flux &a)
//    {
//        return Flux(a.mass*beta, a.momentum_normal*beta,
//                    a.momentum_tangent*beta, a.energy*beta);
//    }
//    friend Flux operator / (const Flux &a, const Set::Scalar beta)
//    {
//        return Flux(a.mass/beta, a.momentum_normal/beta,
//                    a.momentum_tangent/beta, a.energy/beta);
//    }
//
//};

struct Flux {
    std::vector<double> mass;
    Set::Scalar momentum_normal;
    Set::Scalar momentum_tangent;
    Set::Scalar energy;
    Flux() : mass(std::vector<double>(0)), momentum_normal(0.0), momentum_tangent(0.0), energy(0.0)  {}
    Flux(std::vector<double> a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy) :
        mass(a_mass), momentum_normal(a_momentum_normal), momentum_tangent(a_momentum_tangent), energy(a_energy) {}

    Flux(const State &in) :
        mass(in.rho),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
        energy(in.E) {}

    friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
    {
        os << "nspecies=" << flux.mass.size() << ", ";
        for (size_t n=0; n<flux.mass.size(); ++n) {
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
        for (size_t n=0; n<mass.size(); ++n) { mass[n] += a.mass[n]; }
        momentum_normal += a.momentum_normal; momentum_tangent += a.momentum_tangent; 
        energy += a.energy;
    }
    void operator -= (const Flux &a)
    {
        for (size_t n=0; n<mass.size(); ++n) { mass[n] -= a.mass[n]; }
        momentum_normal -= a.momentum_normal; momentum_tangent -= a.momentum_tangent; 
        energy -= a.energy;
    }
    void operator *= (const Set::Scalar alpha)
    {
        for (size_t n=0; n<mass.size(); ++n) { mass[n] *= alpha; }
        momentum_normal *= alpha;
        momentum_tangent *= alpha;
        energy *= alpha;
    }
    void operator /= (const Set::Scalar alpha)
    {
        for (size_t n=0; n<mass.size(); ++n) { mass[n] /= alpha; }
        momentum_normal /= alpha;
        momentum_tangent /= alpha;
        energy /= alpha;
    }
    friend Flux operator + (const Flux &a, const Flux &b)
    {
        std::vector<double> mass_sum(a.mass.size());
        std::transform(a.mass.begin(), a.mass.end(), b.mass.begin(), mass_sum.begin(), std::plus<double>());
        return Flux(mass_sum, a.momentum_normal+b.momentum_normal,
                    a.momentum_tangent+b.momentum_tangent, a.energy+b.energy);
    }
    friend Flux operator - (const Flux &a, const Flux &b)
    {
        std::vector<double> mass_diff(a.mass.size());
        std::transform(a.mass.begin(), a.mass.end(), b.mass.begin(), mass_diff.begin(), std::minus<double>());
        return Flux(mass_diff, a.momentum_normal-b.momentum_normal,
                    a.momentum_tangent-b.momentum_tangent, a.energy-b.energy);
    }
    friend Flux operator * (const Flux &a, const Set::Scalar beta)
    {
        std::vector<double> mass_prod(a.mass.size());
        std::transform(a.mass.begin(), a.mass.end(), mass_prod.begin(), [beta](double val){ return val*beta; });
        return Flux(mass_prod, a.momentum_normal*beta,
                    a.momentum_tangent*beta, a.energy*beta);
    }
    friend Flux operator * (const Set::Scalar beta,const Flux &a)
    {
        std::vector<double> mass_prod(a.mass.size());
        std::transform(a.mass.begin(), a.mass.end(), mass_prod.begin(), [beta](double val){ return val*beta; });
        return Flux(mass_prod, a.momentum_normal*beta,
                    a.momentum_tangent*beta, a.energy*beta);
    }
    friend Flux operator / (const Flux &a, const Set::Scalar beta)
    {
        std::vector<double> mass_quot(a.mass.size());
        std::transform(a.mass.begin(), a.mass.end(), mass_quot.begin(), [beta](double val){ return val/beta; });
        return Flux(mass_quot, a.momentum_normal/beta,
                    a.momentum_tangent/beta, a.energy/beta);
    }

};


class Riemann 
{
public:
//    virtual Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small) = 0;
    virtual Flux Solve(State lo, State hi, Set::Scalar p_ref, Set::Scalar small, std::vector<double> &species_mw, std::vector<double> &species_cp) = 0;
    virtual ~Riemann() = default;
};


}
}
}




#endif
