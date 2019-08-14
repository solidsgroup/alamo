#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string>
#include "Util/Util.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"
#include "Solver/Local/CG.H"
using namespace Model::Solid::CrystalPlastic;

void savefile(float s,float es,const char* loc)
{
	std::fstream myfile;
  	myfile.open (loc, std::ios::out | std::ios::app); 

	myfile << s << "," << es << std::endl;
 	myfile.close();
}
void savefile(std::array<double,12> g, double time,const char* loc)
{
	std::fstream myfile;
  	myfile.open (loc, std::ios::out | std::ios::app); 

	for(int i = 0; i < 12; i++)
	{
		if(i < 11)
		{
			myfile << g[i] << ",";
		}
		else 
		{
			myfile << g[i] << "," << time << std::endl;
		}
	}
	//myfile << std::endl;
 	myfile.close();
}
void deletefile(const char* loc)
{
	if( remove(loc) != 0 )
    perror( "Error deleting file" );
 	else
    puts( "successfully deleted" );
}
/*
Values for the elastic tensor, C:

	K=Bulk Modulus, G=Shear Modulus
	Copper -> K = 123x10^9 Pa (N/m^2) G = 45x10^9 Pa
	L = K - 2/3*G = 9.3x10^10 Pa
	mu = 2/3(K-L) = 2x10^10 Pa
	C11 = L+2*mu, C12 = L, C44 = mu
	C11 = 1.33x10^11 Pa -> 1.33e2 GPa
	C12 = 9.3x10^10 Pa  ->  9.3e GPa
	C44 = 2.0x10^10 Pa  ->  2.0e GPa
 */
int main (int argc, char* argv[])
{
	auto src = static_cast<const char*>("/home/icrman/Python/data_n1.dat");
	auto src2 = static_cast<const char*>("/home/icrman/Python/gamma.dat");
	deletefile(src); deletefile(src2);
	Util::Initialize(argc,argv);
	CrystalPlastic cp(1.33e2, 9.3e1, 2.0e1); //0.4, 0.1, 0.01
	std::array<double,12> gamma;
	//cp.Randomize();
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();
	Set::Matrix sig = Set::Matrix::Zero();
	Set::iMatrix mask = Set::iMatrix::Zero();
	mask(0,0) = 1;
	
	double dt = 1e-8;
	double c = 10.0;
	double T = 0.15/c;
	cp.Setdt(dt);
	int counter = 0;
	for(double t = 0.0; t <= T; t += dt)
	{
		Set::Matrix esp = cp.GetEsp();
		Set::Matrix temp;
		es(0,0) = c*t;
		Eigen::Matrix<amrex::Real,AMREX_SPACEDIM-1,1> s = cp.relax(es, es(0,0));
		es(1,1) = s(0);
		es(2,2) = s(1);

		//temp = es ;
		//es = Solver::Local::CG(cp.DDW(temp),sig,temp,mask,false);

		cp.update(es,sigma,dt);

		if( counter % 1000 == 0)
		{
			for(int i = 0; i < 12; i++)
			{
				gamma[i] = cp.getGamma(i);
			}
			Util::Message(INFO,"es() = ", es);
			savefile( (float)sigma(0,0), (float)es(0,0), src); 
			savefile(gamma,t,src2);
		}
		counter++;
	}
	esp = cp.GetEsp();
	Util::Message(INFO,"esp = ", esp.trace());
	Util::Message(INFO,"sig = ", sigma);
	Util::Message(INFO,"Done");

Util::Finalize();
}