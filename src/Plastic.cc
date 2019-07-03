#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>

#include "Util/Util.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"
using namespace Model::Solid::CrystalPlastic;

void savefile(std::vector<float> s,std::vector<float> es)
{
	std::fstream myfile;
  	myfile.open ("/home/icrman/Python/data.dat",std::ios::out);
	for(int i=0; i < s.size(); i++)
	{
		myfile << s.at(i) << "," << es.at(i) << std::endl;
	}
  	
 	myfile.close();
}

int main (int argc, char* argv[])
{
	CrystalPlastic cp;
	Util::Initialize(argc,argv);
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix s = Set::Matrix::Zero();
	std::vector<float> sig;
	std::vector<float> ees;

	int fin = 70000000;
	for(int i = 0; i < fin; i++)
	{
		es(0,0) = i*1e-8;
	
		cp.SetEs(es);	
		//cp.UpdateSigma();
		for(int j = 0; j <= 1; j++)
		{
			cp.UpdateSigma();
			cp.AdvanceEsp();
		}
		if( i % 10000 == 0)
		{
			sigma = cp.GetSigma();
			sig.emplace_back( (float)sigma(0,0) );
			ees.emplace_back( (float)es(0,0) );
			
			//Util::Message(INFO,"sigma = ", sig.back() );
			
		}
		
	}
	s = cp.GetEsp();
	savefile(sig,ees);
	Util::Message(INFO,"\n",s);
	
	Util::Message(INFO,"Done");
	Util::Finalize();
} 
