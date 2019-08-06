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
	double dt = 1e-6;
	double T = dt;
	cp.Setdt(dt);
	int counter = 0;
	for(double t = 0.0; t <= T; t += dt)
	{
		double c = 1;
		es(0,0) = 5+t*c;
		Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> s = cp.reflux(es, es(0,0));
		es(1,1) = s(0);
		es(2,2) = s(1);
		//es(1,1) = -0.417*t*c;
		//es(2,2) = -0.417*t*c;

		cp.update(es,sigma,dt);

		if( counter % 10000 == 0)
		{
			Util::Message(INFO,"esp() = ", s);
			sig.emplace_back( (float)sigma(0,0) );
			ees.emplace_back( (float)es(0,0) );

		}
		counter++;
	}
	savefile(sig,ees);
	esp = cp.GetEsp();
	Util::Message(INFO,"esp = ", esp.trace());
	Util::Message(INFO,"sig = ", sigma);
/* 
	CrystalPlastic cp(.4,.1,.01);
	Util::Initialize(argc,argv);

	Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> x = Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1>::Zero();
	Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> g;
	x(1) = -1; 
	g = cp.DFP(x, 1e-5 ,0.1 ,0.8 ,1e-5);
	Util::Message(INFO,g);
	Util::Message(INFO,"Done");
*/
	Util::Finalize();
} 
