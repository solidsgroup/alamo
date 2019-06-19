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

	Util::Message(INFO,"Working...");
	Set::Matrix es = Set::Matrix::Zero();
	es(0,0) = 1;
	
	cp.SetEs(es);	
	cp.UpdateSigma();
	cp.AdvanceEsp();
	
	Util::Message(INFO,"Done");

	Util::Finalize();
} 
