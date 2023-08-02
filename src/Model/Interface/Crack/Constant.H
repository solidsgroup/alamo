//
// This class implements an isotropic model of crack.
// This class inherits from the base class :code:`Crack.H`.
// See `this link <https://doi.org/10.1016/j.cma.2021.114011>`_
// for more information on the constant crack model.
//
//
// Class methods:
// 
// #. :code:`Constant()`: 
//    Basic constructor. Does nothing, and leaves all values initiated as NAN.
// #. :code:`Set::Scalar Gc (Set::Scalar /*theta*/)`
//    Returns the value of fracture energy :math:`G(c)`
// #. :code:`Set::Scalar DGc (Set::Scalar /*theta*/)`
//    Returns zero for the constant isotropic crack model
// #. :code:`Set::Scalar DDGc (Set::Scalar /*theta*/)`
//    Returns zero for the constant isotropic crack model
// #. :code:`Set::Scalar Zeta(Set::Scalar /*theta*/)`
//    Returns the crack length scale :math:`\xi`
// #. :code:`Set::Scalar Mobility(Set::Scalar /*theta*/)`
//    Returns the crack mobility :math:`M`
// #. :code:`Set::Scalar DrivingForceThreshold(Set::Scalar /*theta*/)`
//    Returns the driving force threshold.
// #. :code:`static void Parse(Constant & value, IO::ParmParse & pp)`
//    Parses the input file to get model parameters
//
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
        pp.query("G_c",value._Gc);  // crack fracture energy :math:`G_c`
        pp.query("zeta",value.zeta);    // crack length scale :math:`\xi`
        pp.query("mobility",value.mobility);    // crack mobility :math:`M`
        pp.query("threshold", value.threshold); // crack driving force threshold
        
        std::string gtype = "";
        std::string wtype = "";
        pp.query("gtype", gtype);   // choice of degradation function :math:`g(c)`
        pp.query("wtype", wtype);   // choice of interpolation function :math:`w(c)`

        Set::Scalar ductile_exponent = 1.;
        if(gtype == "squarepexp"|| gtype == "cubicm") pp.query("exponent",ductile_exponent);

        std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
        g_map["square"] = Model::Interface::Crack::Crack::GType::GSQUARE;
        g_map["multiwell"] = Model::Interface::Crack::Crack::GType::GMULTIWELL;
        g_map["4c3"] = Model::Interface::Crack::Crack::GType::GPHI4C3;
        g_map["3c2"] = Model::Interface::Crack::Crack::GType::GPHI3C2;
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
