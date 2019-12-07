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
	Crack() {};
	//virtual amrex::Real W(amrex::Real theta) = 0;
	//virtual amrex::Real DW(amrex::Real theta) = 0;
	//virtual amrex::Real DDW(amrex::Real theta) = 0;
	virtual Set::Scalar w_phi(Set::Scalar c)
	{
		switch (m_type)
		{
			case PhiSq:	return 1. - c*c;
			case PhiMultiWell: return (1.-c)*(1.-c)*c*c;
			case Phi4c3: return 1. - 4.*c*c*c + 3.*c*c*c*c;
			default: return 1.-c*c;
		}
	}
	virtual Set::Scalar g_phi(Set::Scalar c)
	{
		switch(m_type)
		{
			case PhiSq:	return c*c;
			case PhiMultiWell: return c*c;
			case Phi4c3: return 4.*c*c*c - 3.*c*c*c*c;
			default: return c*c;
		}
	}
	virtual Set::Scalar Dw_phi(Set::Scalar c)
	{
		switch (m_type)
		{
			case PhiSq:	return -2.*c;
			case PhiMultiWell: return 2.0*c*(1.-c)*(1.-2.*c);
			case Phi4c3: return -12.*c*c*(1.-c);;
			default: return -2.*c;
		}
	}
	virtual Set::Scalar Dg_phi(Set::Scalar c)
	{
		switch (m_type)
		{
			case PhiSq:	return 2.*c;
			case PhiMultiWell: return 2.0*c;
			case Phi4c3: return 12.*c*c*(1.-c);;
			default: return 2.*c;
		}
	}
	//virtual Set::Scalar g_phi(Set::Scalar c) = 0;
	//virtual Set::Scalar Dw_phi(Set::Scalar c) = 0;
	//virtual Set::Scalar Dg_phi(Set::Scalar c) = 0;
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

protected:
	static constexpr amrex::Real pi = 3.14159265359;
	PhiType m_type = PhiType::PhiSq;
};
}
}
}

#endif