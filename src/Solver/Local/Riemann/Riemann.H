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
    std::vector<double>& rho = NAN;
    Set::Scalar M_normal = NAN;
    Set::Scalar M_tangent = NAN;
    Set::Scalar E = NAN;
    // Construtor for convenience
    State() { rho = std::vector<double>; M_normal = 0.0, M_tangent = 0.0; E = 0.0; nspecies = 0; }
    State(std::vector<double>& a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E, int a_nspecies )
        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E), nspecies(a_nspecies) {}
    State(Set::Patch<const Set::Scalar> density_mf,
          Set::Patch<const Set::Scalar> momentum_mf,
          Set::Patch<const Set::Scalar> energy_mf,
          int i, int j, int k, int direction)
    {
        nspecies = density_mf.nComp();
        std::cout << "Riemann nspecies: " << nspecies << std::endl;
        std::cout << "+++++++++++++++++++++++++++++++++" << std::endl;
        for (int n=0; n<nspecies; ++n) { rho[n] = rho_mf(i,j,k,n); }
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

    friend std::ostream &operator<<(std::ostream &os, const State_Extended &state) 
    {
        os << "rho=" << state.rho << ", "; 
        os << "Mn=" << state.M_normal << ", ";
        os << "Mt=" << state.M_tangent << ", ";
        os << "E=" << state.E << ", ";
        os << "nspecies=" << state.nspecies << ", ";
        for (int ii=0; ii<state.nspecies; ++ii) {
            os << "Y[" << ii << "]=" << state.rhoY[ii]/state.rho << ", ";
        }
        // do stuf
        return os;
    }
    void operator += (const State_Extended &a)
    {
        rho += a.rho; M_normal += a.M_normal; M_tangent += a.M_tangent; E += a.E;
        for (int ii=0; ii<nspecies; ++ii) { rhoY[ii] += a.rhoY[ii]; }
    };
    void operator -= (const State_Extended &a)
    {
        rho -= a.rho; M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E; 
        for (int ii=0; ii<nspecies; ++ii) { rhoY[ii] -= a.rhoY[ii]; }
    };
    void operator *= (const Set::Scalar alpha)
    {
        rho *= alpha; M_normal *= alpha; M_tangent *= alpha; E *= alpha; 
        for (int ii=0; ii<nspecies; ++ii) { rhoY[ii] *= alpha; }
    };
    void operator /= (const Set::Scalar alpha)
    {
        rho /= alpha; M_normal /= alpha; M_tangent /= alpha; E /= alpha; 
        for (int ii=0; ii<nspecies; ++ii) { rhoY[ii] /= alpha; }
    };
    friend State_Extended operator + (const State_Extended &a, const State_Extended &b)
    {
        std::vector<double> rhoYsum(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), b.rhoY.begin(), rhoYsum.begin(), std::plus<double>());
        return State_Extended(a.rho+b.rho, a.M_normal+b.M_normal, a.M_tangent+b.M_tangent, a.E+b.E, a.nspecies, rhoYsum);
    }
    friend State_Extended operator - (const State_Extended &a, const State_Extended &b)
    {
        std::vector<double> rhoYdiff(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), b.rhoY.begin(), rhoYdiff.begin(), std::minus<double>());
        return State_Extended(a.rho-b.rho, a.M_normal-b.M_normal, a.M_tangent-b.M_tangent, a.E-b.E, a.nspecies, rhoYdiff);
    }
    friend State_Extended operator * (const Set::Scalar alpha, const State_Extended &b)
    {
        std::vector<double> rhoYprod(b.rhoY.size());
        std::transform(b.rhoY.begin(), b.rhoY.end(), rhoYprod.begin(), [alpha](double val){ return val*alpha; });
        return State_Extended(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha, b.nspecies, rhoYprod);
    }
    friend State_Extended operator * (const State_Extended &b, const Set::Scalar alpha)
    {
        std::vector<double> rhoYprod(b.rhoY.size());
        std::transform(b.rhoY.begin(), b.rhoY.end(), rhoYprod.begin(), [alpha](double val){ return val*alpha; });
        return State_Extended(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha, b.nspecies, rhoYprod);
    }
    friend State_Extended operator / (const State_Extended &b, const Set::Scalar alpha)
    {
        std::vector<double> rhoYquot(b.rhoY.size());
        std::transform(b.rhoY.begin(), b.rhoY.end(), rhoYquot.begin(), [alpha](double val){ return val/alpha; });
        return State_Extended(b.rho/alpha, b.M_normal/alpha, b.M_tangent/alpha, b.E/alpha, b.nspecies, rhoYquot);
    }
    friend State_Extended operator - (const State_Extended &a, const State &b)
    {
        return State_Extended(a.rho-b.rho, a.M_normal-b.M_normal, a.M_tangent-b.M_tangent, a.E-b.E, a.nspecies, a.rhoY);
    }
};

struct Flux {
    Set::Scalar mass;
    Set::Scalar momentum_normal;
    Set::Scalar momentum_tangent;
    Set::Scalar energy;
    Flux() : mass(0.0), momentum_normal(0.0), momentum_tangent(0.0), energy(0.0)  {}
    Flux(Set::Scalar a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy) :
        mass(a_mass), momentum_normal(a_momentum_normal), momentum_tangent(a_momentum_tangent), energy(a_energy) {}

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

struct Flux_Extended {
    Set::Scalar mass;
    Set::Scalar momentum_normal;
    Set::Scalar momentum_tangent;
    Set::Scalar energy;
    int nspecies;
    std::vector<double>& rhoY;
    Flux_Extended() : mass(0.0), momentum_normal(0.0), momentum_tangent(0.0), energy(0.0), rhoY(std::vector<double>(0))  {}
    Flux_Extended(Set::Scalar a_mass, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy, int a_nspecies, std::vector<double>& a_rhoY) :
        mass(a_mass), momentum_normal(a_momentum_normal),
        momentum_tangent(a_momentum_tangent), energy(a_energy), nspecies(a_nspecies), rhoY(a_rhoY) {}

    Flux_Extended(const State_Extended &in) :
        mass(in.rho),
        momentum_normal(in.M_normal),
        momentum_tangent(in.M_tangent),
        energy(in.E),
        nspecies(in.nspecies),
        rhoY(in.rhoY) {}

    friend std::ostream &operator<<(std::ostream &os, const Flux_Extended &flux) 
    {
        os << "mass=" << flux.mass << ", "; 
        os << "Mn=" << flux.momentum_normal << ", ";
        os << "Mt=" << flux.momentum_tangent << ", ";
        os << "E=" << flux.energy << ", ";
        os << "nspecies=" << flux.nspecies << ", ";
        for (int ii=0; ii<flux.nspecies; ++ii) {
            os << "rhoY[" << ii << "]=" << flux.rhoY[ii] << ", ";
        }
        // do stuff
        return os;
    }
    void operator += (const Flux_Extended &a)
    {
        mass += a.mass;
        momentum_normal += a.momentum_normal; momentum_tangent += a.momentum_tangent; 
        energy += a.energy;
        for (int ii=0; ii<a.nspecies; ++ii) { rhoY[ii] += a.rhoY[ii]; }
    }
    void operator -= (const Flux_Extended &a)
    {
        mass -= a.mass;
        momentum_normal -= a.momentum_normal; momentum_tangent -= a.momentum_tangent; 
        energy -= a.energy;
        for (int ii=0; ii<a.nspecies; ++ii) { rhoY[ii] -= a.rhoY[ii]; }
    }
    void operator *= (const Set::Scalar alpha)
    {
        mass *= alpha;
        momentum_normal *= alpha;
        momentum_tangent *= alpha;
        energy *= alpha;
        for (int ii=0; ii<nspecies; ++ii) { rhoY[ii] *= alpha; }
    }
    void operator /= (const Set::Scalar alpha)
    {
        mass /= alpha;
        momentum_normal /= alpha;
        momentum_tangent /= alpha;
        energy /= alpha;
        for (int ii=0; ii<nspecies; ++ii) { rhoY[ii] /= alpha; }
    }
    friend Flux_Extended operator + (const Flux_Extended &a, const Flux_Extended &b)
    {
        std::vector<double> rhoYsum(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), b.rhoY.begin(), rhoYsum.begin(), std::plus<double>());
        return Flux_Extended(a.mass+b.mass, a.momentum_normal+b.momentum_normal,
                    a.momentum_tangent+b.momentum_tangent, a.energy+b.energy,
                    a.nspecies, rhoYsum);
    }
    friend Flux_Extended operator - (const Flux_Extended &a, const Flux_Extended &b)
    {
        std::vector<double> rhoYdiff(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), b.rhoY.begin(), rhoYdiff.begin(), std::minus<double>());
        return Flux_Extended(a.mass-b.mass, a.momentum_normal-b.momentum_normal,
                    a.momentum_tangent-b.momentum_tangent, a.energy-b.energy,
                    a.nspecies, rhoYdiff);
    }
    friend Flux_Extended operator * (const Flux_Extended &a, const Set::Scalar beta)
    {
        std::vector<double> rhoYprod(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), rhoYprod.begin(), [beta](double val){ return val*beta; });
        return Flux_Extended(a.mass*beta, a.momentum_normal*beta,
                    a.momentum_tangent*beta, a.energy*beta,
                    a.nspecies, rhoYprod);
    }
    friend Flux_Extended operator * (const Set::Scalar beta,const Flux_Extended &a)
    {
        std::vector<double> rhoYprod(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), rhoYprod.begin(), [beta](double val){ return val*beta; });
        return Flux_Extended(a.mass*beta, a.momentum_normal*beta,
                    a.momentum_tangent*beta, a.energy*beta,
                    a.nspecies, rhoYprod);
    }
    friend Flux_Extended operator / (const Flux_Extended &a, const Set::Scalar beta)
    {
        std::vector<double> rhoYquot(a.rhoY.size());
        std::transform(a.rhoY.begin(), a.rhoY.end(), rhoYquot.begin(), [beta](double val){ return val/beta; });
        return Flux_Extended(a.mass/beta, a.momentum_normal/beta,
                    a.momentum_tangent/beta, a.energy/beta,
                    a.nspecies, rhoYquot);
    }

};


class Riemann 
{
public:
    //virtual Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small) = 0;
    virtual Flux_Extended Solve(State_Extended lo, State_Extended hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small) = 0;
    virtual ~Riemann() = default;
};


}
}
}




#endif
