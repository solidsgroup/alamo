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
	Crack() {};
	//virtual amrex::Real W(amrex::Real theta) = 0;
	//virtual amrex::Real DW(amrex::Real theta) = 0;
	//virtual amrex::Real DDW(amrex::Real theta) = 0;
	virtual Set::Scalar w_phi(Set::Scalar c) = 0;
	virtual Set::Scalar Dw_phi(Set::Scalar c) = 0;
	virtual Set::Scalar Epc(Set::Scalar theta) = 0;
	virtual Set::Scalar kappa(Set::Scalar theta) = 0;

	void ExportToFile(std::string filename, amrex::Real dTheta)
	{
		std::ofstream outFile;
		outFile.open(filename);
    
		for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
			{
				outFile << theta << " " << W(theta) << std::endl;   
			}
		outFile.close();
    
	}

protected:
	static constexpr amrex::Real pi = 3.14159265359;
};
}
}
}

#endif