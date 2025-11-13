#ifndef MODEL_GAS_EOS_TPG_H_
#define MODEL_GAS_EOS_TPG_H_

#include "Model/Gas/Gas.H"

namespace Model {
namespace Gas {
namespace EOS {

class TPG : public virtual Gas {
public:
    TPG() {};
    virtual ~TPG() = default;

public:
    double ComputeT(double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol=1e-12) const override {
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
            double cv = cv_mass(T, X, i, j, k);
            h = enthalpy_mass(T, X, i, j, k);
            e = h - R(X, i, j, k) * T;
            T = T_old - (e - e_int)/cv;
            convergedT = ( std::fabs(T-T_old)/T < rtol );
        }
        if (!convergedT) {
            Util::Message(INFO, "Temperature didn't converge after ",counter," iterations. T_old: ",T_old," T: ",T);
            Util::Abort(INFO);
        }
        return T;
    }
    double ComputeP(double density, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Pressure, Pa
        double P = density * R(X, i, j, k) * T;
        return P;
    }
    double ComputeE(double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const override {
        // Energy, J/m^3
        double P = density * R(X, i, j, k) * T;
        double rhoE = P / (gamma(T, X, i, j, k) - 1.0);
        double E = rhoE + 0.5*(momentumx*momentumx + momentumy*momentumy)/density;
        return E;
    }

} // class TPG

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
