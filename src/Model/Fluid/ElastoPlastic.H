#ifndef MODEL_ELASTOPLASTIC_H_                                                                                   
#define MODEL_ELASTOPLASTIC_H_                                                                                   

#include <AMReX.H>                                                                                       
#include <AMReX_REAL.H>                                                                                  
#include <eigen3/Eigen/Core>                                                                             
#include "Set/Set.H" 
#include <cmath>                                                                                         
#include <limits>  // For NaN
    
namespace Model {
namespace ElastoPlastic {

class ElastoPlastic {
public:
    ElastoPlastic() {}

    virtual ~ElastoPlastic() {}

    template <typename RealType>
    inline RealType ComputePressureFromDensityAndInternalEnergy(
        RealType rho, RealType e) const {
        if (rho <= 0 || e < 0) {
            Util::Message(INFO, "Invalid values in ComputePressureFromDensityAndInternalEnergy: " + std::to_string(rho));
            return 0.0;
        }
/*
        //Mie-Gruneisen
        Set::Scalar Rhoref = 8930.0;
	Set::Scalar Co = 3940.0;
	Set::Scalar Sref = 1.49;
	Set::Scalar gamma0 = 2.0;

	Set::Scalar gamma = gamma0*Rhoref/rho;
	Set::Scalar d1 = rho*e*gamma;
	Set::Scalar phi = 1.0 - Rhoref/rho;
	Set::Scalar d2 = 1.0 - Sref*phi;
	d2 = std::pow(d2,2);
	d2 = Rhoref*Co*Co*phi/d2;
	Set::Scalar d3 = 1.0 - 0.5*gamma*(rho/Rhoref - 1.0);

        return d1 + d2*d3;
*/
	Set::Scalar Co = 4600.0;
	Set::Scalar Rhoref = 7900.0;
	Set::Scalar gamma0 = 2.17;

	Set::Scalar P = Co*Co*(rho - Rhoref) + (gamma0 - 1.0)*rho*e;

	if (P <= -3.26e8) {
		P = -3.26e8;
	}
	return P;
    }


    template <typename RealType>
    inline RealType ComputeWaveSpeed(
        RealType rho, RealType p, RealType e) const {
        if (rho <= 0 || e < 0) {
            Util::Message(INFO, "Invalid values in ComputeWaveSpeed: " + std::to_string(e));
            return 0.0;
        }
/*
	Set::Scalar Rhoref = 8930.0;
        Set::Scalar Co = 3940.0;
        Set::Scalar Sref = 1.49;
        Set::Scalar gamma0 = 2.0;

        Set::Scalar gamma = gamma0*Rhoref/rho;

	Set::Scalar phi = 1.0 - Rhoref/rho;
	Set::Scalar phi_prime = Rhoref/(rho*rho);
	Set::Scalar d1 = 1.0 - Sref*phi;
	Set::Scalar d3 = std::pow(d1, 3.0);
	Set::Scalar Co2 = std::pow(Co, 2.0);
	Set::Scalar N1 = Rhoref*Co2*phi_prime*(d1 + 2.0*Sref*phi);
	Set::Scalar N2 = 2.0*rho*phi_prime*(d1 + Sref*phi) + phi*d1;
	N2 = -0.5*gamma*Co2*phi*N2;
	Set::Scalar fprime = (N1 + N2)/d3;
	Set::Scalar transverse2 = (4.0/3.0)*(4.5e10/rho);

        return std::sqrt(gamma*e + fprime + transverse2 + gamma*p/rho);
*/
	Set::Scalar Co = 4600.0;
        Set::Scalar Rhoref = 7900.0;
        Set::Scalar gamma0 = 2.17;

	return std::sqrt((gamma0*p + Co*Co*Rhoref)/rho);

    }

    template <typename RealType>
    inline RealType ComputeInternalEnergyFromDensityAndPressure(
        RealType rho, RealType p) const {
        if (rho <= 0) {
            Util::Message(INFO, "Invalid values in ComputeInternalEnergyFromDensityAndPressure: " + std::to_string(rho));
            return 0.0;
        }
/*	Set::Scalar gamma0 = 2.0;
	Set::Scalar Rhoref = 8930.0;
        Set::Scalar Co = 3940.0;
        Set::Scalar Sref = 1.49;

	Set::Scalar gamma = gamma0*Rhoref/rho;
	Set::Scalar d1 = p/(gamma*rho);
	Set::Scalar phi = 1.0 - (Rhoref/rho);
	Set::Scalar d2 = Rhoref*Co*Co*phi;
	d2 = d2/((1.0 - Sref*phi)*(1.0 - Sref*phi));
	Set::Scalar d3 = 1.0 - 0.5*gamma*rho*(1.0/Rhoref - 1.0/rho);

        return d1 - d2*d3/(gamma*rho);
*/
	Set::Scalar Co = 4600.0;
        Set::Scalar Rhoref = 7900.0;
        Set::Scalar gamma0 = 2.17;

	return (p - Co*Co*(rho-Rhoref))/(rho*(gamma0 - 1.0));
    }


    virtual bool ContainsNan() {
        Util::Message(INFO, "Checking for NaN values in ElastoPlastic: " + std::to_string(std::numeric_limits<double>::quiet_NaN()));
        return false;  // Default behavior assumes no NaNs unless overridden.
    }
};

}  // namespace ElastoPlastic
}  // namespace Model

#endif  // MODEL_ELASTOPLASTIC_H_

