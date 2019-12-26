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
    Constant()
    {
        amrex::ParmParse pp("crack.constant");
        pp.query("G_c",Gc);
		pp.query("zeta",zeta);
		pp.query("mobility",mobility);
		pp.query("mult_Gc",mult_Gc);
		pp.query("mult_Lap", mult_Lap);
        std::string phitype = "";
        pp.query("phitype", phitype);

        Set::Scalar ductile_exponent = 1.;
        if(phitype == "squarepexp") pp.query("exponent",ductile_exponent);

        std::map<std::string,Model::Interface::Crack::Crack::PhiType>  phi_map;
	    phi_map["square"] = Model::Interface::Crack::Crack::PhiType::PhiSq;
	    phi_map["multiwell"] = Model::Interface::Crack::Crack::PhiType::PhiMultiWell;
	    phi_map["4c3"] = Model::Interface::Crack::Crack::PhiType::Phi4c3;
        phi_map["squarep"] = Model::Interface::Crack::Crack::PhiType::PhiSqP;
	    phi_map["squarepexp"] = Model::Interface::Crack::Crack::PhiType::PhiSqPM;

        SetPhiType(phi_map[phitype]);
        SetDuctileExponent(ductile_exponent);
    };
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

private:
    Set::Scalar Gc = 1.0, zeta = 1.e-2, mobility = 1.e-2;
    Set::Scalar mult_Gc = 1.0, mult_Lap = 1.0;
};
}
}
}
#endif
