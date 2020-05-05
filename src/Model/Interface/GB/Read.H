#ifndef MODEL_INTERFACE_GB_READ_H
#define MODEL_INTERFACE_GB_READ_H

#include <iostream>
#include <fstream>
#include <vector>

#include "AMReX.H"
#include "GB.H"

#include "Numeric/Interpolator/Linear.H"

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GB
{
/// Reads the data from a file and computes energies and its derivates
///
class Read : public GB
{
public:
	/// \brief Read in data
	///
	/// Reads the data from a file and abort if it is not possible to open the file or if the range of thetas do not give a range between 0 and 2pi. It also computes the derivatives of the energy and stores them in vectors that will be used in W, DW, and DDW
	///
	/// \f[ \int_0^1x^2dx = \frac{1}{3} \f]
	///
	Read() {}
	Read(std::string filename)
	{
		Define(filename);
	}

	void Define(std::string filename)
	{
		std::ifstream input;
		input.open(filename);
		std::string line;
		std::vector<Set::Scalar> theta, thetasmall, w, dw, ddw;
		while(std::getline(input,line))
		{
			std::vector<std::string> dat = Util::String::Split(line);
			theta.push_back(std::stof(dat[0]));
			w.push_back(std::stof(dat[1]));
			Util::Message(INFO,theta[theta.size()-1]," ",w[theta.size()-1]);
		}
		for (unsigned int i = 1; i < theta.size()-1; i++)
		{
			thetasmall.push_back(theta[i]);
			dw.push_back((w[i+1] - w[i-1]) / (theta[i+1] - theta[i-1]));
			ddw.push_back((w[i+1] - 2.0*w[i] + w[i-1]) / ((theta[i+1]-theta[i]) * (theta[i]-theta[i-1])));
		}
		m_w.define(w,theta);
		m_dw.define(dw,thetasmall);
		m_ddw.define(ddw,thetasmall);

		for (Set::Scalar t = -0.001; t < 2*Set::Constant::Pi+1.0; t+=0.001)
		{
			Set::Scalar w = m_w(t), dw = m_dw(t), ddw = m_ddw(t);
			if (std::isnan(w) || std::isnan(dw) || std::isnan(ddw) ||
				std::isinf(w) || std::isinf(dw) || std::isinf(ddw)) 
				Util::Abort(INFO,"Error in GB Read: t=",t," w=",w," dw=",dw," ddw=",ddw);
		}
	};
	Set::Scalar W(amrex::Real theta)
	{
		return m_w(theta);
	};
	Set::Scalar DW(amrex::Real theta)
	{
		return m_dw(theta);
	};
	Set::Scalar DDW(amrex::Real theta)
	{
		return m_ddw(theta);
	};

private:
	Numeric::Interpolator::Linear<Set::Scalar> m_w, m_dw, m_ddw;
	  
public:
	static void Parse(Read & value, amrex::ParmParse & pp)
	{
		std::string filename;
		pp.query("filename",filename);
		value.Define(filename);
	}
	  
};
}
}
}
#endif
