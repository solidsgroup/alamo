
#include <AMReX.H>
#include "MyTest/MyTest.H"


int main (int argc, char* argv[])
{
  std::cout << std::endl << "HERE" << std::endl << std::endl;

  amrex::Initialize(argc, argv);

  {
    MyTest mytest;
    
    mytest.solve();
    mytest.writePlotfile();
  }

    amrex::Finalize();
}
