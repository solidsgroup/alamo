#ifndef SOLVER_LOCAL_RIEMANN_H
#define SOLVER_LOCAL_RIEMANN_H


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
    Set::Scalar eta = NAN;
    // Construtor for convenience
    State() { rho = 0.0; M_normal = 0.0, M_tangent = 0.0; E = 0.0; eta = 0.0;}
    State(Set::Scalar a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E, Set::Scalar a_eta)
        : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E), eta(a_eta) {}
    friend std::ostream &operator<<(std::ostream &os, const State &state) 
    {
        os << "rho=" << state.rho << ", "; 
        os << "Mn=" << state.M_normal << ", ";
        os << "Mt=" << state.M_tangent << ", ";
        os << "E=" << state.E << ", ";
        os << "eta=" << state.eta;
        // do stuf
        return os;
    }
    void operator += (const State &a)
    {
        rho += a.rho; M_normal += a.M_normal; M_tangent += a.M_tangent; E += a.E; eta += a.eta;
    };
    void operator -= (const State &a)
    {
        rho -= a.rho; M_normal -= a.M_normal; M_tangent -= a.M_tangent; E -= a.E; eta -= a.eta;
    };
    void operator *= (const Set::Scalar alpha)
    {
        rho *= alpha; M_normal *= alpha; M_tangent *= alpha; E *= alpha; eta *= alpha;
    };
    void operator /= (const Set::Scalar alpha)
    {
        rho /= alpha; M_normal /= alpha; M_tangent /= alpha; E /= alpha; eta /= alpha;
    };
    friend State operator + (const State &a, const State &b)
    {
        return State(a.rho+b.rho, a.M_normal+b.M_normal, a.M_tangent+b.M_tangent, a.E+b.E, a.eta+b.eta);
    }
    friend State operator - (const State &a, const State &b)
    {
        return State(a.rho-b.rho, a.M_normal-b.M_normal, a.M_tangent-b.M_tangent, a.E-b.E, a.eta-b.eta);
    }
    friend State operator * (const Set::Scalar alpha, const State &b)
    {
        return State(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha, b.eta*alpha);
    }
    friend State operator * (const State &b, const Set::Scalar alpha)
    {
        return State(b.rho*alpha, b.M_normal*alpha, b.M_tangent*alpha, b.E*alpha, b.eta*alpha);
    }
    friend State operator / (const State &b, const Set::Scalar alpha)
    {
        return State(b.rho/alpha, b.M_normal/alpha, b.M_tangent/alpha, b.E/alpha, b.eta/alpha);
    }
};

struct Flux {
    Set::Scalar mass;
    Set::Scalar momentum_normal;
    Set::Scalar momentum_tangent;
    Set::Scalar energy;
    friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
    {
        os << "mass=" << flux.mass << ", "; 
        os << "Mn=" << flux.momentum_normal << ", ";
        os << "Mt=" << flux.momentum_tangent << ", ";
        os << "E=" << flux.energy << ", ";
        // do stuff
        return os;
    }
};



}
}
}




#endif
