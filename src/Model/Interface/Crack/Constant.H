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
        pp.query("G_c",value._Gc);
        pp.query("zeta",value.zeta);
        pp.query("mobility",value.mobility);
        pp.query("threshold", value.threshold);
        
        std::string gtype = "";
        std::string wtype = "";
        pp.query("gtype", gtype);
        pp.query("wtype", wtype);

        Set::Scalar ductile_exponent = 1.;
        if(gtype == "squarepexp"|| gtype == "cubicm") pp.query("exponent",ductile_exponent);

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
