#ifndef MODEL_INTERFACE_GB_READ_H
#define MODEL_INTERFACE_GB_READ_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"

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

	void Define(std::string fileName)
	{
		//Define
		int i;
		std::string line;
		amrex::Real data1, data2; 
		std::ifstream inFile;  

	    //Open file
	    inFile.open(fileName);       
	
	    //Error Check
	    if(inFile.fail()){
	      Util::Abort (INFO, "UNABLE TO OPEN THE FILE");
	    }

	    //Read the file and fill out the vectors thetas and sigmasW with the read values
	    while(getline(inFile, line))
		{ 
	    	int pos = line.find_first_of(' ');
	    	std::string second=line.substr(pos+1);
	    	std::string first=line.substr(0, pos);
	
	    	data1=stod(first);
	    	data2=stod(second);
	
	    	thetas.push_back(data1);
	    	sigmasW.push_back(data2);
	    }

	    size= thetas.size();

	    //Close file
	    inFile.close();

	    //Check that the file we have read contains the whole range of thetas (from 0 to 2pi)
	    if((thetas[size-1]+(thetas[1]-thetas[0]))<2*pi){
	      Util::Abort(INFO, "FILE DID NOT CONTAIN THE WHOLE RANGES OF THETAS (FROM 0 TO 2PI)");
	    }
	
	    //Create vector sigmas DW    
	    for(i=0; i<size; i++)
		{
			if(i==(size-1))
			{
	       		sigmasDW.push_back((sigmasW[0]-sigmasW[i])/(2*pi + thetas[0]-thetas[i]));
	      	}
	      	else
			{
	  			sigmasDW.push_back((sigmasW[i+1]-sigmasW[i])/(thetas[i+1]-thetas[i]));
	    	}
	  	}

	    //Create vector sigmas DDW
	    for(i=0; i<size; i++)
		{
			if(i==(size-1))
			{
				sigmasDDW.push_back((sigmasDW[0]-sigmasDW[i])/(2*pi + thetas[0]-thetas[i]));
			}
	      	else
			{
				sigmasDDW.push_back((sigmasDW[i+1]-sigmasDW[i])/(thetas[i+1]-thetas[i]));
			}
	    }
	    //Print Thetas Sigmas Sigmas' and Sigmas''
	};
	amrex::Real W(amrex::Real theta)
	{
	    // put in code so that 0 <= theta < 2*pi
	    if(theta>=2*pi)
	    {
			int div=theta/(2*pi);
			theta=theta-(div*2*pi);
		}
	    else if(theta<0)
		{
			int div=fabs(theta)/(2*pi);
			theta=theta+((div+1)*2*pi);
	    }
	
	    //Find out the value of theta inside our vector thetas
	
	    for(int i=0; i<size; i++)
		{
			if(i==(size-1))
			{
	  			if(theta >= thetas[i])
				{
	    			//Interpolate to get the value of sigma
	    			return ((theta-thetas[i])*(sigmasW[0]-sigmasW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasW[i];
	  			}
	    	}
	    	else
			{
	  			if(theta >= thetas[i] && theta < thetas[i+1])
				{
	    			//Interpolate to get the value of sigma	  
	    			return ((theta-thetas[i])*(sigmasW[i+1]-sigmasW[i])/(thetas[i+1]-thetas[i]))+sigmasW[i];
	  			}
	    	} 
	    }
	    return 0;
	};
	amrex::Real DW(amrex::Real theta)
	{
		// put in code so that 0 <= theta < 2*pi
	    if(theta>=2*pi)
		{
			int div=theta/(2*pi);
			theta=theta-(div*2*pi);
	    }
	    else if(theta<0)
	    {
			int div=fabs(theta)/(2*pi);
			theta=theta+((div+1)*2*pi);
		}
	
	    //Find out the value of theta inside our vector thetas
	
	    for(int i=0; i<size; i++)
		{
	    	if(i==(size-1))
			{
				if(theta >= thetas[i])
				{
					//Interpolate to get the value of sigma
					return ((theta-thetas[i])*(sigmasDW[0]-sigmasDW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasDW[i];
	  			}
	      	}
	      	else
			{
	  			if(theta >= thetas[i] && theta < thetas[i+1])
				{
	    			//Interpolate to get the value of sigma	  
	    			return ((theta-thetas[i])*(sigmasDW[i+1]-sigmasDW[i])/(thetas[i+1]-thetas[i]))+sigmasDW[i];
	  			}
	    	}
	    }
	    return 0;
	};
	amrex::Real DDW(amrex::Real theta)
	{
		// put in code so that 0 <= theta < 2*pi
		if(theta>=2*pi)
	    {
			int div=theta/(2*pi);
			theta=theta-(div*2*pi);
	    }
	    else if(theta<0)
	    {
			int div=fabs(theta)/(2*pi);
			theta=theta+((div+1)*2*pi);
	    }
	
	    //Find out the value of theta inside our vector thetas
	
	    for(int i=0; i<size; i++)
		{
			if(i==(size-1))
			{
				if(theta >= thetas[i])
				{
					//Interpolate to get the value of sigma
					return ((theta-thetas[i])*(sigmasDDW[0]-sigmasDDW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasDDW[i];
				}
			}
	      	else
			{
	  			if(theta >= thetas[i] && theta < thetas[i+1])
				{
	    			//Interpolate to get the value of sigma	  
	    			return ((theta-thetas[i])*(sigmasDDW[i+1]-sigmasDDW[i])/(thetas[i+1]-thetas[i]))+sigmasDDW[i];
	  			}
	      	} 
	    }
	    return 0;
	};

private:
  	std::vector<amrex::Real> thetas, sigmasW, sigmasDW, sigmasDDW;
  	int size;
	  
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
