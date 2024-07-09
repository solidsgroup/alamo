#ifndef MODEL_INTERFACE_CRACK_SIN_H
#define MODEL_INTERFACE_CRACK_SIN_H

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
class Sin : public Crack
{
public:
    Sin(){};

    Set::Scalar Gc(Set::Scalar theta) {return Gc0 + 0.5*Gc1*(1.0 - cos(4.0*(theta-theta0))); }
    Set::Scalar DGc(Set::Scalar theta) { return 2.0*Gc1*sin(4.0*(theta-theta0)); }
    Set::Scalar DDGc(Set::Scalar theta) { return 8.0*Gc1*cos(4.0*(theta-theta0)); }
    Set::Scalar Zeta(Set::Scalar /*theta*/) { return zeta; }
    Set::Scalar Mobility(Set::Scalar /*theta*/) { return mobility; } 
    Set::Scalar DrivingForceThreshold(Set::Scalar /*theta*/) { return threshold; }

private:
    Set::Scalar zeta = 1.e-2, mobility = 1.e-2;
    Set::Scalar Gc0 = 1.0;
    Set::Scalar Gc1 = 0.49;
    Set::Scalar theta0 = 3.14;
    Set::Scalar threshold = 0.0;

public:
    static void Parse(Sin & value, IO::ParmParse & pp)
    {
        pp_query("Gc0",value.Gc0); // Min Gc (fracture energy) value
        pp_query("Gc1",value.Gc1); // Max Gc (fracture energy)value
        pp_query("theta0",value.theta0); // Angle offset
        pp_query("zeta",value.zeta); // Regularization
        pp_query("mobility",value.mobility); // Crack mobiilty
        pp_query("threshold",value.threshold); // Threshold for kinetics

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
