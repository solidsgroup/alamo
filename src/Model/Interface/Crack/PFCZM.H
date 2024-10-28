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
    Set::Scalar Gc (Set::Scalar /*theta*/) { return _Gc / c_alpha; }
    Set::Scalar DGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar DDGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar Zeta(Set::Scalar /*theta*/) { return zeta; }
    Set::Scalar Mobility (Set::Scalar /*theta*/) {return mobility;}
    Set::Scalar DrivingForceThreshold (Set::Scalar /*theta*/) {return threshold;}

private:
    Set::Scalar _Gc = 1.0e3; 
    Set::Scalar c_alpha = pi; // scaling factor.
    Set::Scalar zeta = 1.e-5;
    Set::Scalar mobility = 1.e-2;
    Set::Scalar threshold = 0.0;
    Set::Scalar czm_a0 = 1.0;
    Set::Scalar czm_order = 1.0;

public:
    static void Parse(PFCZM & value, IO::ParmParse & pp)
    {
        pp_query("czm_order", value.czm_order); //Order 1, 1.5, or 2 (Wu 2024, JMPS)
        pp_query("G_c",value._Gc); // Fracture energy in N/m
        pp_query("zeta",value.zeta); // Lengthscale regularization in m

        Set::Scalar fracture_strength = 1.e6;
        pp_query("f_t", fracture_strength);   // fracture strength beyond which softening starts in Pa
        
        Set::Scalar youngs = 2.e11;
        pp_query("E", youngs);
        
        Set::Scalar irwing_length = youngs * value._Gc / (fracture_strength * fracture_strength);
        value.czm_a0 = 2.0 * irwing_length / (pi * value.zeta);
        Util::Message(INFO, "irwing_length = ", irwing_length, ", a0 = ", value.czm_a0);

        pp_query("mobility",value.mobility); // Mobility (speed)
        pp_query("threshold", value.threshold); // Threshold

        
        std::string gtype = "";
        std::string g2type = "";
        std::string wtype = "";
        pp_query_validate("gtype", gtype, {"wu_linear"}); // degradation function. For now, we only have linear softening
        pp_query_validate("wtype", wtype, {"wu"}); // geometric function. no choice except for the optimal one for now.
        pp_query_validate("g2type", g2type, {"wu_linear"}); // dissipation function. For now, we only have linear softening

        std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
        g_map["wu_linear"] = Model::Interface::Crack::Crack::GType::GWULINEAR;
        
        std::map<std::string,Model::Interface::Crack::Crack::WType>  w_map;
        w_map["wu"] = Model::Interface::Crack::Crack::WType::WU;

        std::map<std::string,Model::Interface::Crack::Crack::G2Type>  g2_map;
        g2_map["wu_linear"] = Model::Interface::Crack::Crack::G2Type::G2WULINEAR;
        
        value.SetGType(g_map[gtype]);
        value.SetWType(w_map[wtype]);
        value.SetG2Type(g2_map[g2type]);

        if (w_map[wtype] == Model::Interface::Crack::Crack::WType::WU) value.c_alpha = pi;

        value.SetPFCZMConstants(value.czm_order, value.czm_a0);
    }
};
}
}
}
#endif
