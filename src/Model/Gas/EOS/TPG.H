#ifndef MODEL_GAS_EOS_TPG_H_
#define MODEL_GAS_EOS_TPG_H_

#include <cmath>
#include <memory>
#include "Model/Gas/EOS/EOS.H"
#include "Util/Util.H"

namespace Model {
namespace Gas {
namespace EOS {

class TPG : public EOS {
public:
    TPG() {} ;
    //TPG(std::shared_ptr<GasData> data_in, std::shared_ptr<Gas> gas_in) : EOS(data_in, gas_in) {};
    ~TPG() override = default;

public:
    double ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol=1e-12) const override {
        /*
        // Temperature, K
        double vx = momentumx/density;
        double vy = momentumy/density;
        double e_int = E/density - 0.5*(vx*vx + vy*vy);

        // Newton-Raphson Method
        double T_old = Tguess;
        double T = Tguess;
        bool convergedT = false;
        double h=0.0, e=0.0;
        double counter = 0;
        while (!convergedT && counter <= 100) {
            ++counter;
            T_old = T;
            double cv = gas->cv_mass(T, X, i, j, k);
            h = gas->enthalpy_mass(T, X, i, j, k);
            e = h - gas->R(X, i, j, k) * T;
            T = T_old - (e - e_int)/cv;
            convergedT = ( std::fabs(T-T_old)/T < rtol );
        }
        if (!convergedT) {
            Util::Abort(INFO, "Temperature didn't converge after ",counter," iterations. T_old: ",T_old," T: ",T);
        }
        return T;
        */
    }
    double ComputeT(double pressure, double density,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        /*
        // Temperature, K
        double T = pressure / density / gas->R(X, i, j, k);
        return T;
        */
    }
    double ComputeP(double density, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        /*
        // Pressure, Pa
        double P = density * gas->R(X, i, j, k) * T;
        return P;
        */
    }
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        /*
        // Energy, J/m^3
        double P = density * gas->R(X, i, j, k) * T;
        double rhoE = P / (gas->gamma(T, X, i, j, k) - 1.0);
        double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
        */
    }

}; // class TPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
