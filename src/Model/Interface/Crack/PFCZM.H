#ifndef MODEL_INTERFACE_CRACK_PFCZM_H
#define MODEL_INTERFACE_CRACK_PFCZM_H

#include <iostream>
#include <fstream>
#include <cmath>

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
    enum FSType {WANG2023, MC};

    PFCZM() = default;
    virtual ~PFCZM() = default;

    // In Wu 2024 formulation, Gc is scaled by a constant C_alpha.
    // $C_\alpha = 4 \int_0^1 \sqrt{\alpha(\beta)} d\beta$
    // where $\alpha$ is the geometric function defined in the code as w_phi
    // For the optimal function chosen in Wu 2024, $C_\alpha = \pi$ 
    Set::Scalar Gc (Set::Scalar /*ratio*/) { return _Gc / _c_alpha; }
    Set::Scalar c_alpha ()      { return _c_alpha; }
    Set::Scalar l_w()           { return _l_w; }

    // this hould check if we need to do mixed mode loading.
    bool mixed_mode()           { return _mixed_mode; }
    FSType failure_surface()    { return _failure_surface; }

    // added functions for mixed mode - Wang 2023 surface
    Set::Scalar GcI_bar ()      { return _GcI_bar / _c_alpha; }
    Set::Scalar GcII_bar ()     { return _GcII_bar / _c_alpha; }
    Set::Scalar chi()           { return _chi; }
    Set::Scalar beta_bar()      { return _beta_bar; }
    Set::Scalar sig_t()         { return _sig_t; }
    Set::Scalar tau_s()         { return _tau_s; }

    // added functions for mixed mode - mohr columb surface
    Set::Scalar cohesion()      { return _cohesion; }
    Set::Scalar friction()      { return _friction; }

    // For now, we are storing things like youngs modulus stuff here.
    // ideally, they should be moved to Material -> Solid -> Affine.
    Set::Scalar DGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar DDGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar Zeta(Set::Scalar /*theta*/) { return zeta; }
    Set::Scalar Mobility (Set::Scalar /*theta*/) {return mobility;}
    Set::Scalar DrivingForceThreshold (Set::Scalar /*theta*/) {return threshold;}

private:
    Set::Scalar _Gc = NAN; 
    Set::Scalar _c_alpha = pi; // scaling factor.
    Set::Scalar _l_w = 1.e-5; //irwing length
    Set::Scalar zeta = NAN;
    Set::Scalar mobility = NAN;
    Set::Scalar threshold = NAN;

    Set::Scalar czm_a0 = 1.0;
    Set::Scalar czm_order = NAN;
    
    bool _mixed_mode = false;
    FSType _failure_surface = FSType::WANG2023;

    Set::Scalar _GcI_bar = 1.e3;
    Set::Scalar _GcII_bar = 1.e3;

    // Parameters specific to WANG 2023 failure surface
    Set::Scalar _chi = NAN;
    Set::Scalar k1 = 1.0;
    Set::Scalar _beta_bar = -1.0;
    Set::Scalar _sig_t = 1.0;
    Set::Scalar _tau_s = 1.0;

    // Parameters specific to MC failure surface
    Set::Scalar _cohesion = NAN;
    Set::Scalar _friction = NAN;

public:
    static void Parse(PFCZM & value, IO::ParmParse & pp)
    {
        pp.query_default("czm_order", value.czm_order, 1.0); //Order 1, 1.5, or 2 (Wu 2024, JMPS)
        pp.query_default("G_c",value._Gc,1.0e3); // Fracture energy in N/m
        pp.query_default("zeta",value.zeta,1.e-5); // Lengthscale regularization in m

        Set::Scalar fracture_strength = NAN;
        
        Set::Scalar youngs = NAN;
        Set::Scalar nu = NAN;
        Set::Scalar tensile_strength = 2.e8;
        Set::Scalar shear_strength = 1.41e8;
        Set::Scalar compressive_strength = NAN;

        // Toggle to determine mixed mode
        pp.query_default("mixed_mode", value._mixed_mode, false);

        if (value._mixed_mode)
        {
            // Young's modulus
            pp.query_default("E", youngs, 2.e11);
            
            // Poisson's ratio $\nu$
            pp.query_default("nu", nu, 0.3);

            // shear modulus
            Set::Scalar mu = 0.5 * youngs / (1. + nu);
            
            std::string fail_surf_type = "";
            // type of failure surface
            pp.query_validate("failure_surface", fail_surf_type, {"wang2023", "mohr", "mc", "columb"}); // degradation function. For now, we only have linear softening
            std::map<std::string,FSType>  fs_map;
            fs_map["wang2023"]   = FSType::WANG2023;
            fs_map["mohr"]       = FSType::MC;
            fs_map["mc"]         = FSType::MC;
            fs_map["columb"]     = FSType::MC;
            value._failure_surface = fs_map[fail_surf_type];

            if (value._failure_surface == FSType::WANG2023)
            {
                // Fracture strength
                pp.query_default("fracture_strength", fracture_strength, 1.e6);
                
                // $\chi$ value
                pp.query_default("chi", value._chi, 1.0);
                
                // Compressive strength
                pp.query_default("compressive_strength", compressive_strength,4.e8);
                value.k1 = ( value._chi < 1/std::sqrt(2) ) ? ( 2 * value._chi * std::sqrt(1.0 - value._chi*value._chi) ) : 1.0;

                tensile_strength = fracture_strength / value.k1;
                shear_strength = value._chi * tensile_strength;
                value._beta_bar = (1.0 / (value._chi * value._chi)) - (compressive_strength * compressive_strength / (4.0 * value._chi * value._chi * shear_strength * shear_strength));

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
            else if (value._failure_surface == FSType::MC)
            {
                // cohesion coefficient for MC FS type
                pp.query_default("cohesion", value._cohesion, 1.0);
                // friction coefficient for MC FS type
                pp.query_default("friction", value._friction, 0.5);
                
                Set::Scalar phi = std::atan(value._friction);
                fracture_strength = 2.0 * value._cohesion * std::cos(phi) / (1.0 + std::sin(phi));
            }

        }
        else
        {
            // Young's modulus
            pp.query_default("E", youngs, 2.e11);
            // fracture strength beyond which softening starts [Pa]
            pp.query_default("fracture_strength", fracture_strength,1.e6);   
        }
        
        Set::Scalar irwing_length = youngs * value._Gc / (fracture_strength * fracture_strength);
        value.czm_a0 = 4.0 * irwing_length / (pi * value.zeta);
        Util::Message(INFO, "irwing_length = ", irwing_length, ", a0 = ", value.czm_a0);
        value._l_w = irwing_length;

        if (value.zeta > 0.33 * irwing_length) Util::Warning(INFO, "Zeta value is greater than irwing length. Consider reducing it.");

        pp.query_default("mobility",value.mobility,1.0e-2); // Mobility (speed)
        pp.query_default("threshold", value.threshold,0.0); // Threshold

        std::string gtype = "";
        std::string wtype = "";
        pp.query_validate("gtype", gtype, {"wu_linear","square","4c3"}); // degradation function. For now, we only have linear softening
        pp.query_validate("wtype", wtype, {"wu","square","4c3"}); // geometric function. no choice except for the optimal one for now.
        
        std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
        g_map["wu_linear"] = Model::Interface::Crack::Crack::GType::GWULINEAR;
        g_map["4c3"] = Model::Interface::Crack::Crack::GType::GPHI4C3;
        g_map["square"] = Model::Interface::Crack::Crack::GType::GSQUARE;
        
        std::map<std::string,Model::Interface::Crack::Crack::WType>  w_map;
        w_map["wu"] = Model::Interface::Crack::Crack::WType::WU;
        w_map["4c3"] = Model::Interface::Crack::Crack::WType::WPHI4C3;
        w_map["square"] = Model::Interface::Crack::Crack::WType::WSQUARE;

        value.SetGType(g_map[gtype]);
        value.SetWType(w_map[wtype]);

        if (w_map[wtype] == Model::Interface::Crack::Crack::WType::WU) value._c_alpha = pi;
        else if (w_map[wtype] == Model::Interface::Crack::Crack::WType::WSQUARE) value._c_alpha = 4.0/3.0;
        else if (w_map[wtype] == Model::Interface::Crack::Crack::WType::WPHI4C3) value._c_alpha = 12.0/5.0;

        value.SetPFCZMConstants(value.czm_order, value.czm_a0);
    }
};
}
}
}
#endif
