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
    PFCZM(){};
    // In Wu 2024 formulation, Gc is scaled by a constant C_alpha.
    // $C_\alpha = 4 \int_0^1 \sqrt{\alpha(\beta)} d\beta$
    // where $\alpha$ is the geometric function defined in the code as w_phi
    // For the optimal function chosen in Wu 2024, $C_\alpha = \pi$ 
    Set::Scalar Gc (Set::Scalar ratio) { return _Gc / _c_alpha; }
    Set::Scalar c_alpha () { return _c_alpha; }

    // added functions for mixed mode
    Set::Scalar GcI_bar () {return _GcI_bar / _c_alpha; }
    Set::Scalar GcII_bar () {return _GcII_bar / _c_alpha; }
    Set::Scalar chi() { return _chi; }
    Set::Scalar beta_bar() { return _beta_bar; }

    Set::Scalar sig_t() { return _sig_t; }
    Set::Scalar tau_s() { return _tau_s; }

    // For now, we are storing things like youngs modulus stuff here.
    // ideally, they should be moved to Material -> Solid -> Affine.
    
    bool mixed_mode() { return _mixed_mode; }
    
    Set::Scalar DGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar DDGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar Zeta(Set::Scalar /*theta*/) { return zeta; }
    Set::Scalar Mobility (Set::Scalar /*theta*/) {return mobility;}
    Set::Scalar DrivingForceThreshold (Set::Scalar /*theta*/) {return threshold;}

private:
    Set::Scalar _Gc = 1.0e3; 
    Set::Scalar _c_alpha = pi; // scaling factor.
    Set::Scalar zeta = 1.e-5;
    Set::Scalar mobility = 1.e-2;
    Set::Scalar threshold = 0.0;

    Set::Scalar czm_a0 = 1.0;
    Set::Scalar czm_order = 1.0;
    
    bool _mixed_mode = false;

    Set::Scalar _GcI_bar = 1.e3;
    Set::Scalar _GcII_bar = 1.e3;

    Set::Scalar _chi = 1.0;
    Set::Scalar k1 = 1.0;
    Set::Scalar k2 = 1.0;
    Set::Scalar _beta_bar = -1.0;

    Set::Scalar _sig_t = 1.0;
    Set::Scalar _tau_s = 1.0;

public:
    static void Parse(PFCZM & value, IO::ParmParse & pp)
    {
        pp.query("czm_order", value.czm_order); //Order 1, 1.5, or 2 (Wu 2024, JMPS)
        pp.query("G_c",value._Gc); // Fracture energy in N/m
        pp.query("zeta",value.zeta); // Lengthscale regularization in m

        Set::Scalar fracture_strength = 1.e6;
        
        Set::Scalar youngs = 2.e11;
        Set::Scalar nu = 0.3;
        Set::Scalar tensile_strength = 2.e8;
        Set::Scalar shear_strength = 1.41e8;
        Set::Scalar compressive_strength = 4.e8;

        bool mixed_mode = false;
        pp.query("mixed_mode", value._mixed_mode);

        if (value._mixed_mode)
        {
            pp.query("E", youngs);
            pp.query("nu", nu);
            pp.query("fracture_strength", fracture_strength);
            pp.query("chi", value._chi);
            pp.query("compressive_strength", compressive_strength);
            value.k1 = ( value._chi < 1/std::sqrt(2) ) ? ( 2 * value._chi * std::sqrt(1.0 - value._chi*value._chi) ) : 1.0;

            tensile_strength = fracture_strength / value.k1;
            shear_strength = value._chi * tensile_strength;
            value._beta_bar = (1.0 / (value._chi * value._chi)) - (compressive_strength * compressive_strength / (4.0 * value._chi * value._chi * shear_strength * shear_strength));

            

            // pp.query("tensile_strength", tensile_strength);
            // pp.query("shear_strength", shear_strength);
            // pp.query("compressive_strength", compressive_strength);
            // Set::Scalar chi = shear_strength/tensile_strength;

            Set::Scalar mu = 0.5 * youngs / (1. + nu);

            // value._chi = chi;
            // value.k1 = ( chi < 1/std::sqrt(2) ) ? ( 2 * chi * std::sqrt(1.0 - chi*chi) ) : 1.0;
            // value._beta_bar = (1.0 / (chi*chi)) - ((compressive_strength*compressive_strength)/(4.0*chi*chi*shear_strength*shear_strength));
            
            value._GcI_bar = value._Gc/(value.k1*value.k1);
            value._GcII_bar = value._GcI_bar * youngs * value._chi * value._chi / mu ;
            fracture_strength = tensile_strength * value.k1;

            Util::Message(INFO, "Chi = ", value._chi);
            Util::Message(INFO, "mu = ", mu);
            Util::Message(INFO, "k1 = ", value.k1);
            Util::Message(INFO, "beta_bar = ", value._beta_bar);
            Util::Message(INFO, "GcI_bar = ", value._GcI_bar);
            Util::Message(INFO, "GcII_bar = ", value._GcII_bar);
            Util::Message(INFO, "f_t = ", fracture_strength);
            Util::Message(INFO, "f_c = ", compressive_strength);
            Util::Message(INFO, "sig_t = ", tensile_strength);
            Util::Message(INFO, "tau_s = ", shear_strength);

            value._sig_t = tensile_strength;
            value._tau_s = shear_strength;

        }
        else
        {
            pp.query("E", youngs);
            pp.query("fracture_strength", fracture_strength);   // fracture strength beyond which softening starts in Pa
        }
        
        Set::Scalar irwing_length = youngs * value._Gc / (fracture_strength * fracture_strength);
        value.czm_a0 = 4.0 * irwing_length / (pi * value.zeta);
        Util::Message(INFO, "irwing_length = ", irwing_length, ", a0 = ", value.czm_a0);

        if (value.zeta > 0.33 * irwing_length) Util::Warning(INFO, "Zeta value is greater than irwing length. Consider reducing it.");

        pp.query("mobility",value.mobility); // Mobility (speed)
        pp.query("threshold", value.threshold); // Threshold

        std::string gtype = "";
        std::string wtype = "";
        pp.query_validate("gtype", gtype, {"wu_linear", "ai_cubic"}); // degradation function. For now, we only have linear softening
        pp.query_validate("wtype", wtype, {"wu"}); // geometric function. no choice except for the optimal one for now.
        
        std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
        g_map["wu_linear"] = Model::Interface::Crack::Crack::GType::GWULINEAR;
        g_map["ai_cubic"] = Model::Interface::Crack::Crack::GType::GAICUBIC;
        
        std::map<std::string,Model::Interface::Crack::Crack::WType>  w_map;
        w_map["wu"] = Model::Interface::Crack::Crack::WType::WU;

        value.SetGType(g_map[gtype]);
        value.SetWType(w_map[wtype]);

        if (w_map[wtype] == Model::Interface::Crack::Crack::WType::WU) value._c_alpha = pi;

        value.SetPFCZMConstants(value.czm_order, value.czm_a0);
    }
};
}
}
}
#endif
