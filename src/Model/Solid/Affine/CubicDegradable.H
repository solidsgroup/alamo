#ifndef MODEL_SOLID_AFFINE_CUBIC_DEGRADABLE_H_
#define MODEL_SOLID_AFFINE_CUBIC_DEGRADABLE_H_

#include "AMReX.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Affine/Affine.H"
#include "Model/Solid/Affine/Cubic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class CubicDegradable : public Cubic
{
public:

    enum LatticeType {FCC, BCC, HCP};   // standard Bravais lattices. More to be added later

    CubicDegradable() {};
    CubicDegradable(Solid<Set::Sym::MajorMinor> base) : Cubic(base) {};
    CubicDegradable(Cubic base) : Cubic(base) {};

    CubicDegradable(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Set::Scalar a_phi1, Set::Scalar a_Phi, Set::Scalar a_phi2)
    {
        Eigen::Matrix3d m;
        m =  Eigen::AngleAxisd(a_phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(a_Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(a_phi1, Eigen::Vector3d::UnitX());
        Define(a_C11, a_C12, a_C44, m);
    };

    void Define(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44,Eigen::Matrix3d m)
    {
        R = m; C110 = a_C11; C120 = a_C12; C440 = a_C44;
        Cubic::Define(a_C11, a_C12, a_C44, m);
    }

    void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		if(std::isnan(C110)) Util::Abort(INFO, "Nans in C110");
		if(std::isnan(C120)) Util::Abort(INFO, "Nans in C120");
        if(std::isnan(C440)) Util::Abort(INFO, "Nans in C440");

		Cubic::Define( C110*(1.-eta), C120*(1.-eta), C440*(1.-eta), R, F0);
	}
    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
	{
        if(eta.size() == 0) return;
        else if (eta.size() == 1) DegradeModulus(eta[0]);
        else if(eta.size() == 2)
		{
			if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            Cubic::Define( C110*(1.-eta[0]), C120*(1.-eta[1]), C440*(1.-eta[1]), R, F0 );
        }
        else if (eta.size() == 3)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1]) || std::isnan(eta[2])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1]) || std::isnan(eta[2])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0 || eta[2] < 0.0 || eta[2] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            Cubic::Define( C110*(1.-eta[0]), C120*(1.-eta[1]), C440*(1.-eta[2]), R, F0);
        }
        else
            return;
    }
    
    void UpdateF0(Set::Matrix a_F0){F0 = a_F0; Cubic::UpdateF0(a_F0);}

public:
    Set::Scalar C110, C120, C440;
    Set::Matrix F0;
    Eigen::Matrix3d R;

public:
    static CubicDegradable Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }

    static CubicDegradable Random(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44)
    {
        CubicDegradable ret;
        Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
	    Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
	    Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
        Eigen::Matrix3d m;
        m = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());

        ret.Define(a_C11,a_C12,a_C44,m);
        return ret;
    }

    static void Parse(CubicDegradable & value, IO::ParmParse & pp)
    {
        Set::Scalar C11 = 1.68, C12 = 1.21, C44 = 0.75;
        pp.query("C11",C11);
        pp.query("C12",C12);
        pp.query("C44",C44);

        //if (pp.contains("random"))
        //{
        //    value = CubicCrystalPlastic::Random(C11,C12,C44);
        //    return;
        //}

        Set::Scalar phi1 = 0.0, Phi = 0.0, phi2 = 0.0;
        pp.query("phi1",phi1);
        pp.query("Phi",Phi);
        pp.query("phi2",phi2);

        Eigen::Matrix3d m;
    	m = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());

        value.Define(C11,C12,C44,m);
    }
};

}
}
}
#endif