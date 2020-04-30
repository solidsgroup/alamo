#ifndef MODEL_SOLID_LINEAR_CUBIC_H_
#define MODEL_SOLID_LINEAR_CUBIC_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Cubic : public Solid<Set::Sym::MajorMinor>
{
public:

	Cubic() {};
    Cubic(Solid<Set::Sym::MajorMinor> base) : Solid<Set::Sym::MajorMinor>(base) {};
	virtual ~Cubic() {};

    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
    {
    	Eigen::Matrix3d m;
    	m =     Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
    	Define(C11,C12,C44,m);
    }
    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44, Eigen::Matrix3d R)
    {
    
    	amrex::Real Ctmp[3][3][3][3];
    	ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();

    	for(int i = 0; i < 3; i++) 
    		for(int j = 0; j < 3; j++) 
    			for(int k = 0; k < 3; k++) 
    				for(int l = 0; l < 3; l++)
    				{
    					if(i == j && j == k && k == l)  Ctmp[i][j][k][l] = C11;
    					else if (i==k && j==l) Ctmp[i][j][k][l] = C44;
    					else if (i==j && k==l) Ctmp[i][j][k][l] = C12;
    					else Ctmp[i][j][k][l] = 0.0;
    				}
    	for(int p = 0; p < AMREX_SPACEDIM; p++) 
    		for(int q = 0; q < AMREX_SPACEDIM; q++) 
    			for(int s = 0; s < AMREX_SPACEDIM; s++) 
    				for(int t = 0; t < AMREX_SPACEDIM; t++)
    				{
    					ddw(p,q,s,t) = 0.0;
    					for(int i = 0; i < 3; i++) 
    						for(int j = 0; j < 3; j++) 
    							for(int k = 0; k < 3; k++) 
    								for(int l = 0; l < 3; l++) 
    									ddw(p,q,s,t) += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
    				}
    }
    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return ddw*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return ddw;
    }
    virtual void Print(std::ostream &out) const //override 
    {
        out << ddw;
    }

public:
    static const KinematicVariable kinvar = KinematicVariable::gradu;


    static Cubic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static Cubic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        Cubic ret;
	    Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
	    Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
	    Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
        ret.Define(C11,C12,C44,phi1,Phi,phi2);
        return ret;
    }

    static void Parse(Cubic & value, IO::ParmParse & pp)
    {
        Set::Scalar C11 = 1.68, C12 = 1.21, C44 = 0.75;
        pp.query("C11",C11);
        pp.query("C12",C12);
        pp.query("C44",C44);

        if (pp.contains("random"))
        {
            value = Cubic::Random(C11,C12,C44);
            return;
        }

        Set::Scalar phi1 = 0.0, Phi = 0.0, phi2 = 0.0;
        pp.query("phi1",phi1);
        pp.query("Phi",Phi);
        pp.query("phi2",phi2);
        value.Define(C11,C12,C44,phi1,Phi,phi2);
    }
};


}
}
}

#endif

