#include <iostream>
#include <fstream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <string>
#include "Util/Util.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"
#include "Integrator/EshelbyPlastic.H"
#include "IO/ParmParse.H"
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
	Util::Initialize(argc,argv);

	amrex::ParmParse pp;
	std::string src = "/home/icrman/Python/data_n0.dat";
	std::string src2 = "/home/icrman/Python/gamma.dat";
	pp.query("src",src);
	pp.query("src2",src2);

	//std::ofstream fsrc;  fsrc.open(src);
	//std::ofstream fsrc2; fsrc2.open(src);

	deletefile(src.c_str()); deletefile(src2.c_str());
	CrystalPlastic cp( 10e3 * 2.75e-3, 7.35e3 * 2.75e-3, 3.8e3 * 2.75e-3); //0.4, 0.1, 0.01 C11 = c11*t0
	// 10e3 * 0.75e-3, 7.35e3 * 0.75e-3, 3.8e3 * 0.75e-3
	std::array<double,12> gamma;
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix sigma = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();
	Set::iMatrix mask = Set::iMatrix::Zero();
	mask(0,1) = 1; mask(1,0) = 1;
	//mask(0,0) = 1;
	Set::Scalar dt = 1e-4;
	Set::Scalar c = 1.0;
	Set::Scalar T = 0.8/c;
	cp.Setdt(dt);
	int counter = 0;
		
	for(double t = 0.0; t <= T; t += dt)
	{
		esp = cp.GetEsp();
		Set::Matrix temp;
		//es(0,0) = c*t; 
		es(0,1) = c*t; es(1,0) = c*t;
		temp = (es - esp);
		es = Solver::Local::CG(cp.DDW(es),-sigma,es,mask,false);
		cp.update(es,sigma);

		if( counter % 10 == 0)
		{
			gamma = cp.StressSlipSystem(sigma);
			for(int i = 0; i < 12; i++)
			{
				gamma[i] = cp.getGamma(i);
			}
			//Util::Message(INFO,"es() = ", es);
			//Util::Message(INFO,"esp() = ", esp);
			//Util::Message(INFO,"sig() = ", sigma);
			//Util::Message(INFO,"t=",t," es(0,1)=",es(0,1));
			//fsrc << es(0,0) << " " << sigma(0,0) << std::endl;
			//fsrc1 << es(0,1) << "," << sigma(0,1) << std::endl;
			savefile((float)es(0,1), (float)sigma(0,1), src.c_str()); 
			savefile(gamma,t,src2.c_str());
		}
		counter++;
	}
	Util::Message(INFO,"es() = ", es, "\ntr(esp) = ", esp.trace());
	Util::Message(INFO,"esp() = ", esp);
	Util::Message(INFO,"sig() = ", sigma);

	//fsrc.close();
	Util::Finalize(); return 0;
}