#include <iostream>

#if AMREX_SPACEDIM==2

#include <AMReX.H>
#include "MyTest/MyTest.H"


int main (int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  {
    MyTest mytest;
    
    mytest.solve();
    mytest.writePlotfile();
  }

    amrex::Finalize();
}
#else
int main()
{
    std::cout << "This program works in 2D only" << std::endl;
}
#endif
