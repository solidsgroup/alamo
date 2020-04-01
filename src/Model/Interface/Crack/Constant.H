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
    Set::Scalar Epc(Set::Scalar /*theta*/)
    {
        return mult_Gc*Gc/(4.0*zeta);
    };
    Set::Scalar kappa(Set::Scalar /*theta*/)
    {
        return mult_Lap*2.0*Gc*zeta;
    };
    Set::Scalar GetMobility( Set::Scalar /*theta*/)
    {
        return mobility;
    }
	Set::Scalar DrivingForceThreshold(Set::Scalar /*theta*/)
	{
		return threshold;
	}

private:
    Set::Scalar Gc = 1.0, zeta = 1.e-2, mobility = 1.e-2, threshold = 0.0;
    Set::Scalar mult_Gc = 1.0, mult_Lap = 1.0;

public:
    static void Parse(Constant & value, IO::ParmParse & pp)
    {
        pp.query("G_c",value.Gc);
		pp.query("zeta",value.zeta);
		pp.query("mobility",value.mobility);
		pp.query("mult_Gc",value.mult_Gc);
		pp.query("mult_Lap", value.mult_Lap);
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
