#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "Model/Solid/LinearElastic/CrystalPlastic.H"
using namespace Model::Solid::CrystalPlastic;

int main (int argc, char* argv[])
{
	CrystalPlastic cp;	
	Util::Initialize(argc,argv);
	Set::Matrix s;
	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	Set::Matrix esp = Set::Matrix::Zero();
	int fin = 2/1e-7;
	for(int i = 0; i < 100; i++)
	{
		es(0,0) = 1e-3*i;
	
		cp.SetEs(es);	
		//cp.SetEsp(esp);
		cp.UpdateSigma();
		for(int j = 0; j < 100000; j++)
		{
			esp = cp.AdvanceEsp();
		}
		s = cp.GetEsp();
		Util::Message(INFO,"esp = ", s);
	}
	
	Util::Message(INFO,"Done");

	Util::Finalize();
} 
