#ifndef MODEL_INTERFACE_CRACK_CONSTANT_H
#define MODEL_INTERFACE_CRACK_CONSTANT_H

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
class Constant : public Crack
{
public:
    Constant(){};
    
    Set::Scalar Gc (Set::Scalar /*theta*/) { return _Gc; }
    Set::Scalar DGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar DDGc (Set::Scalar /*theta*/) { return 0.0; }
    Set::Scalar Zeta(Set::Scalar /*theta*/) { return zeta; }
    Set::Scalar Mobility (Set::Scalar /*theta*/) {return mobility;}
    Set::Scalar DrivingForceThreshold (Set::Scalar /*theta*/) {return threshold;}
private:
    Set::Scalar _Gc = 1.0, zeta = 1.e-2, mobility = 1.e-2, threshold = 0.0;

public:
    static void Parse(Constant & value, IO::ParmParse & pp)
    {
        pp_query("G_c",value._Gc); // Fracture energy 
        pp_query("zeta",value.zeta); // Lengthscale regularization
        pp_query("mobility",value.mobility); // Mobility (speed)
        pp_query("threshold", value.threshold); // Threshold
        
        std::string gtype = "";
        std::string wtype = "";
        pp_query("gtype", gtype); // Type of g function to use {square, multiwell, 4c3, squarep, squarepexp, cubicm}
        pp_query("wtype", wtype); // Type o w function to use {square, multiwell, multiwell2, 4c3}

        Set::Scalar ductile_exponent = 1.;
        if(gtype == "squarepexp"|| gtype == "cubicm") pp_query("exponent",ductile_exponent);

        std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
        g_map["square"] = Model::Interface::Crack::Crack::GType::GSQUARE;
        g_map["multiwell"] = Model::Interface::Crack::Crack::GType::GMULTIWELL;
        g_map["4c3"] = Model::Interface::Crack::Crack::GType::GPHI4C3;
        g_map["squarep"] = Model::Interface::Crack::Crack::GType::GSQUAREP;
        g_map["squarepexp"] = Model::Interface::Crack::Crack::GType::GSQUAREPM;
        g_map["cubicm"] = Model::Interface::Crack::Crack::GType::GCUBICM;

        std::map<std::string,Model::Interface::Crack::Crack::WType>  w_map;
        w_map["square"] = Model::Interface::Crack::Crack::WType::WSQUARE;
        w_map["multiwell"] = Model::Interface::Crack::Crack::WType::WMULTIWELL;
        w_map["multiwell2"] = Model::Interface::Crack::Crack::WType::WMULTIWELL2;
        w_map["4c3"] = Model::Interface::Crack::Crack::WType::WPHI4C3;

        value.SetGType(g_map[gtype]);
        value.SetWType(w_map[wtype]);

        value.SetDuctileExponent(ductile_exponent);
    }
};
}
}
}
#endif
