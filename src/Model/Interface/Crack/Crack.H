#ifndef MODEL_INTERFACE_CRACK_H
#define MODEL_INTERFACE_CRACK_H

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <iostream>
#include <fstream>

namespace Model
{
namespace Interface
{
namespace Crack
{
class Crack
{
	public:
	enum PhiType {PhiSq, PhiMultiWell, Phi4c3}; 
	enum DuctilePhiType {PhiSqP, PhiSqPM}

	Crack() {};
	//virtual amrex::Real W(amrex::Real theta) = 0;
	//virtual amrex::Real DW(amrex::Real theta) = 0;
	//virtual amrex::Real DDW(amrex::Real theta) = 0;

	// Methods for brittle fracture
	virtual Set::Scalar w_phi(Set::Scalar c)
	{
		switch (m_type)
		{
			case PhiSq:			return 1. - c*c;
			case PhiMultiWell: 	return (1.-c)*(1.-c)*c*c;
			case Phi4c3: 		return 1. - 4.*c*c*c + 3.*c*c*c*c;
			default: 			return 1.-c*c;
		}
	}
	virtual Set::Scalar g_phi(Set::Scalar c)
	{
		switch(m_type)
		{
			case PhiSq:			return c*c;
			case PhiMultiWell: 	return c*c;
			case Phi4c3: 		return 4.*c*c*c - 3.*c*c*c*c;
			default: 			return c*c;
		}
	}
	virtual Set::Scalar Dw_phi(Set::Scalar c)
	{
		switch (m_type)
		{
			case PhiSq:			return -2.*c;
			case PhiMultiWell: 	return 2.0*c*(1.-c)*(1.-2.*c);
			case Phi4c3: 		return -12.*c*c*(1.-c);;
			default: 			return -2.*c;
		}
	}
	virtual Set::Scalar Dg_phi(Set::Scalar c)
	{
		switch (m_type)
		{
			case PhiSq:			return 2.*c;
			case PhiMultiWell: 	return 2.0*c;
			case Phi4c3: 		return 12.*c*c*(1.-c);;
			default: 			return 2.*c;
		}
	}

	// Methods for ductile fracture
	virtual Set::Scalar w_phi(Set::Scalar c, Set::Scalar p)
	{
		switch (m_type_d)
		{
			case PhiSqP:	return 1. - c*c;
			case PhiSqPM:	return 1. - c*c;
			default: 		return 1. - c*c;
		}
	}
	virtual Set::Scalar Dw_phi(Set::Scalar c, Set::Scalar p)
	{
		switch (m_type_d)
		{
			case PhiSqP:	return -2.*c;
			case PhiSqPM:	return -2.*c;
			default: 		return -2.*c;
		}
	}
	virtual Set::Scalar g_phi(Set::Scalar c, Set::Scalar p)
	{
		switch (m_type_d)
		{
			case PhiSqP:	return std::pow(c,2.*p);
			case PhiSqPM:	return std::pow(c,2.*(std::pow(p,m_d_exponent)));
			default: 		return std::pow(c,2.*p);;
		}
	}
	virtual Set::Scalar Dg_phi(Set::Scalar c, Set::Scalar p)
	{
		switch (m_type_d)
		{
			case PhiSqP:	return 2.*p*std::pow(c,2.*p - 1.);
			case PhiSqPM:	return 2.*std::pow(p,m_d_exponent)*std::pow(c,2.*(std::pow(p,m_d_exponent)-1.));
			default: 		return 2.*p*std::pow(c,2.*p - 1.);
		}
	}

	virtual Set::Scalar Epc(Set::Scalar theta) = 0;
	virtual Set::Scalar kappa(Set::Scalar theta) = 0;

	void ExportToFile(std::string filename, amrex::Real dTheta)
	{
		std::ofstream outFile;
		outFile.open(filename);
    
		for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
			{
				outFile << theta << " " << Epc(theta) << std::endl;   
			}
		outFile.close();
    
	}
	void SetPhiType(const PhiType a_type)
	{
		m_type = a_type;
	}
	void SetDuctilePhiType(const DuctilePhiType a_type)
	{
		m_type_d = a_type;
	}
	void SetDuctileExponent(const Set::Scalar m)
	{
		m_d_exponent = m;
	}

protected:
	static constexpr amrex::Real pi = 3.14159265359;
	PhiType m_type = PhiType::PhiSq;
	DuctilePhiType m_type_d = DuctilePhiType::PhiSqP;
	Set::Scalar m_d_exponent = 1.;
};
}
}
}

#endif