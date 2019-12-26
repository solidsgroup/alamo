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
    Sin()
    {
        amrex::ParmParse pp("crack.sin");
        pp.query("Gc0",Gc0);
        pp.query("Gc1",Gc1);
        pp.query("theta0",theta0);
		pp.query("zeta",zeta);
		pp.query("mobility",mobility);
		pp.query("mult_Gc",mult_Gc);
		pp.query("mult_Lap", mult_Lap);

        std::string phitype = "";
        pp.query("phitype", phitype);

        Set::Scalar ductile_exponent = 1.;
        if(phitype == "squareexp") pp.query("exponent",ductile_exponent);

        std::map<std::string,Model::Interface::Crack::Crack::PhiType>  phi_map;
	    phi_map["square"] = Model::Interface::Crack::Crack::PhiType::PhiSq;
	    phi_map["multiwell"] = Model::Interface::Crack::Crack::PhiType::PhiMultiWell;
	    phi_map["4c3"] = Model::Interface::Crack::Crack::PhiType::Phi4c3;
        phi_map["squarep"] = Model::Interface::Crack::Crack::PhiType::PhiSqP;
	    phi_map["squarepexp"] = Model::Interface::Crack::Crack::PhiType::PhiSqPM;

        SetPhiType(phi_map[phitype]);
        SetDuctileExponent(ductile_exponent);
    };
    Set::Scalar Epc(Set::Scalar theta)
    {
        return mult_Gc*(Gc0 + 0.5*Gc1*(1.0 - cos(4.0*(theta-theta0))))/(4.0*zeta);
    };
    Set::Scalar kappa(Set::Scalar theta)
    {
        return mult_Lap*2.0*(Gc0 + 0.5*Gc1*(1.0 - cos(4.0*(theta-theta0))))*zeta;
    };
    Set::Scalar GetMobility( Set::Scalar /*theta*/)
    {
        return mobility;
    }

private:
    Set::Scalar zeta = 1.e-2, mobility = 1.e-2;
    Set::Scalar mult_Gc = 1.0, mult_Lap = 1.0;
    Set::Scalar Gc0 = 1.0;
    Set::Scalar Gc1 = 0.49;
    Set::Scalar theta0 = 3.14;
};
}
}
}
#endif
