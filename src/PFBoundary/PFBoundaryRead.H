#ifndef PFBoundaryRead_H
#define PFBoundaryRead_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "PFBoundary.H"

#define PI 3.14159265 

class PFBoundaryRead : public PFBoundary
{
public:
  PFBoundaryRead(std::string fileName)
  {
    //Define
    int i;
    std::string line;
    amrex::Real data1, data2;
    std::vector<amrex::Real> thetas, sigmasW, sigmasDW, sigmasDDW;    
    std::ifstream inFile;  

    //Open file
    inFile.open(fileName);       
    
    //Error Check
    if(inFile.fail()){
      std::cout << "\n**UNABLE TO OPEN "+fileName+"**\n" << std::endl;
    }

    //Read the file and fill out the vectors thetas and sigmasW with the read values
    i=0;

    thetas.push_back(0.0);
    sigmasW.push_back(0.0);
    while(getline(inFile, line)){           
      int pos = line.find_first_of(' ');
      std::string second=line.substr(pos+1);
      std::string first=line.substr(0, pos);
      
      data1=stod(first);
      data2=stod(second);
      
      thetas[i]=data1;        
      sigmasW[i]=data2;
      thetas.push_back(0.0);
      sigmasW.push_back(0.0);
      i++;
    }

    thetas.pop_back();
    sigmasW.pop_back();
    
    int size= thetas.size();

    //Close file
    inFile.close();
    
    //Create vector sigmas DW    
    sigmasDW.push_back(0.0);
    for(i=0; i<size; i++){
      if(i==(size-1)){
     	sigmasDW[i]=(sigmasW[0]-sigmasW[i])/(2*pi + thetas[0]-thetas[i]);
     	sigmasDW.push_back(0.0);
      }
      else{
	sigmasDW[i]=(sigmasW[i+1]-sigmasW[i])/(thetas[i+1]-thetas[i]);
      sigmasDW.push_back(0.0);
    }
  }

    sigmasDW.pop_back();

    //Create vector sigmas DDW
    sigmasDDW.push_back(0.0);
    for(i=0; i<size; i++){
      if(i==(size-1)){
	sigmasDDW[i]=(sigmasDW[0]-sigmasDW[i])/(2*pi + thetas[0]-thetas[i]);
	sigmasDDW.push_back(0.0);
      }
      else{
	sigmasDDW[i]=(sigmasDW[i+1]-sigmasDW[i])/(thetas[i+1]-thetas[i]);
	sigmasDDW.push_back(0.0);
      }
    }

    sigmasDDW.pop_back();

    //Print Thetas Sigmas Sigmas' and Sigmas''
    
    std::cout << "Thetas      Sigmas       Sigmas'        Sigmas''" << std::endl;
    for(i=0; i<size; i++){
    std::cout << thetas[i] << "       " << sigmasW[i] << "       " << sigmasDW[i] << "       " << sigmasDDW[i] << std::endl;
    }
  };
  amrex::Real W(amrex::Real theta)
  {
    // put in code so that 0 <= theta < 2*pi
    if(theta>=2*pi)
      {
	theta=theta-2*pi;
      }
    else if(theta<0)
      {
	theta=theta+2*pi;
      }
    
    // use modulo operation on floating point numbers to adjust theta
    
    //Find out the value of theta inside our vector thetas
    
    for(int i=0; i<size; i++){
      if(i==(size-1)){
	if(theta >= thetas[i]){
	  //Interpolate to get the value of sigma
	  return ((theta-thetas[i])*(sigmasW[0]-sigmasW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasW[i];
	}
      }
      else{
	if(theta >= thetas[i] && theta < thetas[i+1]){
	  //Interpolate to get the value of sigma	  
	  return ((theta-thetas[i])*(sigmasW[i+1]-sigmasW[i])/(thetas[i+1]-thetas[i]))+sigmasW[i];
	}
      } 
    }   
  };
  amrex::Real DW(amrex::Real theta)
  {
     // put in code so that 0 <= theta < 2*pi
    if(theta>=2*pi)
      {
	theta=theta-2*pi;
      }
    else if(theta<0)
      {
	theta=theta+2*pi;
      }
    
    // use modulo operation on floating point numbers to adjust theta
    
    //Find out the value of theta inside our vector thetas
    
    for(int i=0; i<size; i++){
      if(i==(size-1)){
	if(theta >= thetas[i]){
	  //Interpolate to get the value of sigma'
	  return ((theta-thetas[i])*(sigmasDW[0]-sigmasDW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasDW[i];
	}
      }
      else{
	if(theta >= thetas[i] && theta < thetas[i+1]){
	  //Interpolate to get the value of sigma'	  
	  return ((theta-thetas[i])*(sigmasDW[i+1]-sigmasDW[i])/(thetas[i+1]-thetas[i]))+sigmasDW[i];
	}
      } 
    }   
  };
  amrex::Real DDW(amrex::Real theta)
  {
     // put in code so that 0 <= theta < 2*pi
    if(theta>=2*pi)
      {
	theta=theta-2*pi;
      }
    else if(theta<0)
      {
	theta=theta+2*pi;
      }
    
    // use modulo operation on floating point numbers to adjust theta
    
   //Find out the value of theta inside our vector thetas
    
    for(int i=0; i<size; i++){
      if(i==(size-1)){
	if(theta >= thetas[i]){
	  //Interpolate to get the value of sigma''
	  return ((theta-thetas[i])*(sigmasDDW[0]-sigmasDDW[i])/(2*pi+thetas[0]-thetas[i]))+sigmasDDW[i];
	}
      }
      else{
	if(theta >= thetas[i] && theta < thetas[i+1]){
	  //Interpolate to get the value of sigma''	  
	  return ((theta-thetas[i])*(sigmasDDW[i+1]-sigmasDDW[i])/(thetas[i+1]-thetas[i]))+sigmasDDW[i];
	}
      } 
    }   
  };

private:
  std::vector<amrex::Real> thetas, sigmasW, sigmasDW, sigmasDDW;
  int size;
};

#endif
