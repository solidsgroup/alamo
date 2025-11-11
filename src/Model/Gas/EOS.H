#ifndef MODEL_GAS_EOS_H
#define MODEL_GAS_EOS_H

#include <AMReX_Array.H>
#include <cmath>
#include "Util/Util.H"
#include "Model/Gas/Species.H"

// EOS functions to compute missing state variables
// Since Model::Gas::Species is templated to return properties as function of T,
// we will just default to the thermally perfect EOS, P=rho*R*T

namespace Model {
namespace Gas {
namespace EOS {

template <typename CompType>
AMREX_GPU_HOST_DEVICE
double ComputeD(const amrex::Array4<double>& rhoY, const int i, const int j, const int k, const int nspecies) {
    double density = 0.0;
    for (int n=0; n<nspecies; ++n) density += rhoY(i,j,k,n);
    return density;
}

template <typename CompType>
AMREX_GPU_HOST_DEVICE
double ComputeT(double density, double momentumx, double momentumy, 
                double E, const CompType& gas) {
    double vx = momentumx/density;
    double vy = momentumy/density;
    double e_int = E/density - 0.5*(vx*vx + vy*vy);

    // Newton-Raphson Method
    double T_old = 1000.0;
    double T = 1000.0;
    bool convergedT = false;
    double rtol = 1e-12;
    double h=0.0, e=0.0;
    double counter = 0;
    while (!convergedT && counter <= 100) {
        ++counter;
        T_old = T;
        double cv = gas.cv_mass(T);
        h = gas.enthalpy_mass(T);
        e = h - gas.R*T;
        T = T_old - (e - e_int)/cv;
        convergedT = ( std::fabs(T-T_old)/T < rtol );
    }
    if (!convergedT) {
        Util::Message(INFO, "Temperature didn't converge after ",counter," iterations. T_old: ",T_old," T: ",T);
        //Util::Abort(INFO);
    }
    return T;
}
    
template <typename CompType>
AMREX_GPU_HOST_DEVICE
double ComputeP(double density, double T, const CompType& gas) {
    return density*gas.R*T;
}

template <typename CompType>
AMREX_GPU_HOST_DEVICE
double ComputeE(double density, double momentumx, double momentumy,
                double T, const CompType& gas) {
    double vx = momentumx/density;
    double vy = momentumy/density;
    double h = gas.enthalpy_mass(T);
    double e = h - gas.R*T;
    double E = density * (e + 0.5*(vx*vx + vy*vy));
    return E;
}

} // namespace EOS
} // namespace Gas
} // namesapce Model

#endif
