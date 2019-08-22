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
	Aluminum -> K = 69x10^9 Pa G = 27x10^9 Pa
	L = K - 2/3*G = 9.3x10^10 Pa
	mu = 2/3(K-L) = 2x10^10 Pa
	C11 = L+2*mu, C12 = L, C44 = mu 
	Copper
	C11 = 1.33x10^11 Pa -> 1.33e2 GPa
	C12 = 9.3x10^10 Pa  ->  9.3e1 GPa
	C44 = 2.0x10^10 Pa  ->  2.0e1 GPa
	Aluminum
	C11 = 7.5x10^10 Pa -> 7.5e1 GPa
	C12 = 2.1x10^10 Pa -> 2.1e1 GPa
	C44 = 1.2x10^10 Pa -> 1.2e1 GPa
 */
int main (int argc, char* argv[])
{
	auto src = static_cast<const char*>("/home/icrman/Python/data_n0.dat");
	auto src2 = static_cast<const char*>("/home/icrman/Python/gamma.dat");
	deletefile(src); deletefile(src2);
	Util::Initialize(argc,argv);
	CrystalPlastic cp(10e3*1.5e-3,7.35e3*1.5e-3,3.8e3*1.5e-3); //0.4, 0.1, 0.01 C11 = c11*t0
	std::array<double,12> gamma;
	//CrystalPlastic cp(0.4,0.1,0.01,.4,.1,.6);
	//cp.Randomize();
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();

	Set::Matrix sig = Set::Matrix::Zero();
	Set::iMatrix mask = Set::iMatrix::Zero();
	mask(0,0) = 1;
	
	double dt = 1e-7;
	double c = 10.0;
	double T = 0.05;
	cp.Setdt(dt);
	int counter = 0;
	for(double t = 0.0; t <= T; t += dt)
	{
		Set::Matrix esp = cp.GetEsp();
		Set::Matrix temp;
		es(0,0) = c*t;
		es = cp.relax(es, es(0,0));
	//	temp = (es - esp);
	//	es = Solver::Local::CG(cp.DDW(es),-cp.DW(temp),es,mask,false);

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