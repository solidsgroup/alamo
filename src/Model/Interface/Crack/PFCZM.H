#ifndef MODEL_INTERFACE_CRACK_PFCZM_H
#define MODEL_INTERFACE_CRACK_PFCZM_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "Crack.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Model
{
namespace Interface
{
namespace Crack
{
class PFCZM : public Crack
{
public:
    PFCZM() = default;
    virtual ~PFCZM() = default;

    // In Wu 2024 formulation, Gc is scaled by a constant C_alpha.
    // $C_\alpha = 4 \int_0^1 \sqrt{\alpha(\beta)} d\beta$
    // where $\alpha$ is the geometric function defined in the code as w_phi
    // For the optimal function chosen in Wu 2024, $C_\alpha = \pi$ 
    Set::Scalar Gc (Set::Scalar _Gc) override { return _Gc / _c_alpha; }
    Set::Scalar c_alpha() { return _c_alpha; }
    Set::Scalar l_w(Set::Scalar _Gc) { return _youngs * _Gc / (_fracture_strength * _fracture_strength); }

    // added functions for mixed mode
    Set::Scalar GcI_bar (Set::Scalar _Gc) {return _Gc / (_k1*_k1*_c_alpha); }
    Set::Scalar GcII_bar (Set::Scalar _GcI_bar) {return _GcI_bar * _youngs * _chi * _chi / (_mu * _c_alpha); }

    Set::Scalar chi() { return _chi; }
    Set::Scalar beta_bar() { return _beta_bar; }

    Set::Scalar sig_t() { return _sig_t; }
    Set::Scalar tau_s() { return _tau_s; }

    // For now, we are storing things like youngs modulus stuff here.
    // ideally, they should be moved to Material -> Solid -> Affine.
    
    bool mixed_mode() { return _mixed_mode; }
    
    Set::Scalar DGc() override { return 0.0; }
    Set::Scalar DDGc() override { return 0.0; }
    Set::Scalar Zeta() override { return _zeta; }
    Set::Scalar Mobility() override { return _mobility; }
    Set::Scalar DrivingForceThreshold() override { return _threshold; }

protected:
    Set::Scalar _m_pf_czm_a0(Set::Scalar _Gc) override { return 4.0 * l_w(_Gc) / (pi * _zeta); }

private:
    Set::Scalar _c_alpha = NAN;
    Set::Scalar _zeta = NAN;
    Set::Scalar _mobility = NAN;
    Set::Scalar _threshold = NAN;

    Set::Scalar _youngs = NAN;
    Set::Scalar _nu = NAN;
    Set::Scalar _tensile_strength = NAN;
    Set::Scalar _shear_strength = NAN;
    Set::Scalar _compressive_strength = NAN;

    Set::Scalar _fracture_strength = NAN;

    Set::Scalar _mu = NAN;
  
    bool _mixed_mode = false;

    Set::Scalar _chi = NAN;
     Set::Scalar _k1 = NAN;
    Set::Scalar _beta_bar = NAN;

    Set::Scalar _sig_t = NAN;
    Set::Scalar _tau_s = NAN;

public:
    static void Parse(PFCZM & value, IO::ParmParse & pp)
    {
        pp.query_default("c_alpha", value._c_alpha, pi);
        pp.query_default("czm_order", value._m_pf_czm_order, 1.0); //Order 1, 1.5, or 2 (Wu 2024, JMPS)
        pp.query_default("zeta", value._zeta, 1e-5); // Lengthscale regularization in m
        
        // Toggle to determine mixed mode
        pp.query("mixed_mode", value._mixed_mode);

        // Young's modulus
        pp.query_default("E", value._youngs, 2e11);
        // Fracture strength
        pp.query_default("fracture_strength", value._fracture_strength, 1e6);

        if (value._mixed_mode)
        {
            // Poisson's ratio $\nu$
            pp.query_default("nu", value._nu, 0.3);
            // $\chi$ value
            pp.query_default("chi", value._chi, 1.0);
            // Compressive strength
            pp.query("compressive_strength", value._compressive_strength, 4e8);

            value._k1 = ( value._chi < 1/std::sqrt(2) ) ? ( 2 * value._chi * std::sqrt(1.0 - value._chi*value._chi) ) : 1.0;

            value._tensile_strength = value._fracture_strength / value._k1;
            value._shear_strength = value._chi * value._tensile_strength;
            value._beta_bar = (1.0 / (value._chi * value._chi)) - (value._compressive_strength * value._compressive_strength / (4.0 * value._chi * value._chi * value._shear_strength * value._shear_strength));
            value._mu = 0.5 * value._youngs / (1. + value._nu);
            
            value._fracture_strength = value._tensile_strength * value._k1;

            Util::Message(INFO, "Chi = ", value._chi);
            Util::Message(INFO, "mu = ", value._mu);
            Util::Message(INFO, "k1 = ", value._k1);
            Util::Message(INFO, "beta_bar = ", value._beta_bar);
            Util::Message(INFO, "f_t = ", value._fracture_strength);
            Util::Message(INFO, "f_c = ", value._compressive_strength);
            Util::Message(INFO, "sig_t = ", value._tensile_strength);
            Util::Message(INFO, "tau_s = ", value._shear_strength);

            value._sig_t = value._tensile_strength;
            value._tau_s = value._shear_strength;

        }

        pp.query_default("mobility", value._mobility, 1e-2); // Mobility (speed)
        pp.query_default("threshold", value._threshold, 0.0); // Threshold

        std::string gtype = "";
        std::string wtype = "";
        pp.query_validate("gtype", gtype, {"wu_linear"}); // degradation function. For now, we only have linear softening
        pp.query_validate("wtype", wtype, {"wu"}); // geometric function. no choice except for the optimal one for now.
        
        std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
        g_map["wu_linear"] = Model::Interface::Crack::Crack::GType::GWULINEAR;
        
        std::map<std::string,Model::Interface::Crack::Crack::WType>  w_map;
        w_map["wu"] = Model::Interface::Crack::Crack::WType::WU;

        value.SetGType(g_map[gtype]);
        value.SetWType(w_map[wtype]);

        if (w_map[wtype] == Model::Interface::Crack::Crack::WType::WU) value._c_alpha = pi;
    }
};
}
}
}
#endif
