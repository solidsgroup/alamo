#include "Viscoelastic.H"
#include <iomanip>

namespace Model
{
namespace Solid
{
namespace Viscoelastic
{
std::ostream&
operator<< (std::ostream& os,
	    const Viscoelastic&    /*b*/)
{
//#if AMREX_SPACEDIM == 2
//	std::array<Set::Matrix,6> gradu, eps;
//	gradu[0] << 1,0, 0,0; eps[0] = 0.5*(gradu[0] + gradu[0].transpose());
//	gradu[1] << 0,0, 0,1; eps[1] = 0.5*(gradu[1] + gradu[1].transpose());
//	gradu[5] << 0,1, 0,0; eps[5] = 0.5*(gradu[5] + gradu[5].transpose());
//	for (int i = -1; i < 7; i++)
//	{
//		if (i==-1) os << "┌                                                               ┐" << std::endl;
//		//else if (i==5) os << "└ ";
//		else  if (i < 6)
//		{
//			os << "│ ";
//			for (int j = 0; j < 6; j++)
//			{
//				if ((i==0 || i==1 || i==5) &&
//				    (j==0 || j==1 || j==5))
//				{
//					Set::Scalar comp = (eps[i].transpose() * b(eps[j])).trace();
//					os << std::setw(10) << (fabs(comp)>1E-10 ? comp : 0);
//				}
//				else
//					os << std::setw(10) << "-";
//			}
//			os << "  │" << std::endl;
//		}
//		else os << "└                                                               ┘" << std::endl;
//
//	}
//#endif
//#if AMREX_SPACEDIM == 3
//	std::array<Set::Matrix,6> gradu, eps;
//	gradu[0] << 1,0,0, 0,0,0, 0,0,0; eps[0] = 0.5*(gradu[0] + gradu[0].transpose());
//	gradu[1] << 0,0,0, 0,1,0, 0,0,0; eps[1] = 0.5*(gradu[1] + gradu[1].transpose());
//	gradu[2] << 0,0,0, 0,0,0, 0,0,1; eps[2] = 0.5*(gradu[2] + gradu[2].transpose());
//	gradu[3] << 0,0,0, 0,0,1, 0,0,0; eps[3] = 0.5*(gradu[3] + gradu[3].transpose());
//	gradu[4] << 0,0,1, 0,0,0, 0,0,0; eps[4] = 0.5*(gradu[4] + gradu[4].transpose());
//	gradu[5] << 0,1,0, 0,0,0, 0,0,0; eps[5] = 0.5*(gradu[5] + gradu[5].transpose());
//	for (int i = -1; i < 7; i++)
//	{
//		if (i==-1) os << "┌                                                               ┐" << std::endl;
//		//else if (i==5) os << "└ ";
//		else  if (i < 6)
//		{
//			os << "│ ";
//			for (int j = 0; j < 6; j++)
//			{
//				Set::Scalar comp = (eps[i].transpose() * b(eps[j])).trace();
//				os << std::setw(10) << (fabs(comp)>1E-10 ? comp : 0);
//			}
//			os << "  │" << std::endl;
//		}
//		else os << "└                                                               ┘" << std::endl;
//
//	}
//#endif
	return os;
}

}
}
}

