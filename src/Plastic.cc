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
  	myfile.open ("/home/icrman/Python/data_n2.dat",std::ios::out);
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
	//CrystalPlastic(1.5,2,4.5);
	//cp.Randomize();
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();;
	
	std::vector<float> sig;
	std::vector<float> ees;
	double T = 10;
	double dt = 1e-6;
	cp.Setdt(dt);
	int counter = 0;
	for(double t = 0.0; t <= T; t+=dt)
	{
		es(0,0) = t;
		sigma = cp.UpdateSigma(es);
		cp.update(es,sigma,dt);

		if( counter % 10000 == 0)
		{
			sig.emplace_back( (float)sigma(0,0) );
			ees.emplace_back( (float)es(0,0) );
			esp = cp.GetEsp();
			//Util::Message(INFO,"esp = ", esp);
			//Util::Message(INFO,"sigma = ", sig.back() );
			
		}
		counter++;
	}
	savefile(sig,ees);
	


	/* 
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
	*/
	
	Util::Message(INFO,"Done");
	Util::Finalize();
} 
