#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string>
#include "Util/Util.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"
using namespace Model::Solid::CrystalPlastic;

void savefile(float s,float es,const char* loc)
{
	std::fstream myfile;
  	myfile.open (loc, std::ios::out | std::ios::app); 

	myfile << s << "," << es << std::endl;
 	myfile.close();
}
void deletefile(const char* loc)
{
	if( remove(loc) != 0 )
    perror( "Error deleting file" );
 	else
    puts( "File successfully deleted" );
}
int main (int argc, char* argv[])
{
	auto src = static_cast<const char*>("/home/icrman/Python/data_n1.dat");
	deletefile(src);
	Util::Initialize(argc,argv);
	CrystalPlastic cp(0.4, 0.1, 0.01);
	//cp.Randomize();
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();
	
	double dt = 1e-7;
	double T = 0.5;
	cp.Setdt(dt);
	int counter = 0;
	for(double t = 0.0; t <= T; t += dt)
	{
		double c = 10;
		es(0,0) = t*c;
		Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> s = cp.relax(es, es(0,0));
		es(1,1) = s(0);
		es(2,2) = s(1);
		//es(1,1) = -0.417*t*c;
		//es(2,2) = -0.417*t*c;

		cp.update(es,sigma,dt);

		if( counter % 10000 == 0)
		{
			Util::Message(INFO,"es() = ", s);
			savefile( (float)sigma(0,0), (float)es(0,0), src); 
		}
		counter++;
	}
	esp = cp.GetEsp();
	Util::Message(INFO,"esp = ", esp.trace());
	Util::Message(INFO,"sig = ", sigma);
	Util::Message(INFO,"Done");

	Util::Finalize();
} 
