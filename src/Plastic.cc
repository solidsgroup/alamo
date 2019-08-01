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
  	myfile.open ("/home/icrman/Python/data_n0.dat",std::ios::out);
	for(unsigned int i=0; i < s.size(); i++)
	{
		myfile << s.at(i) << "," << es.at(i) << std::endl;
	}
  	
 	myfile.close();
}

int main (int argc, char* argv[])
{
	CrystalPlastic cp(.4,.1,.01);
	Util::Initialize(argc,argv);
	//cp.Randomize();
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();
	
	std::vector<float> sig;
	std::vector<float> ees;
	double T = 5;
	double dt = 1e-6;
	cp.Setdt(dt);
	int counter = 0;
	for(double t = 0.0; t <= T; t += dt)
	{
		double c = 1;
		double s = cp.removeStress(es,-t/2);
		es(0,0) = t*c;
		es(1,1) = s;
		es(2,2) = s;
		//es(1,1) = -0.5*t*c;
		//es(2,2) = -0.5*t*c;

		cp.update(es,sigma,dt);

		if( counter % 10000 == 0)
		{
			Util::Message(INFO,"esp(1,1) = ", s);
			sig.emplace_back( (float)sigma(0,0) );
			ees.emplace_back( (float)es(0,0) );
			//esp = cp.GetEsp();
			//Util::Message(INFO,"esp = ", esp);
		}
		counter++;
	}
	savefile(sig,ees);
	
	Util::Message(INFO,"Done");
	Util::Finalize();
} 
