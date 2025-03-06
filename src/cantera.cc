#ifdef ALAMO_CANTERA
#include "cantera/core.h"
#include "cantera/base/AnyMap.h"
#include <cantera/zeroD/IdealGasReactor.h>
#include <cantera/zeroD/ReactorNet.h>
#include <cantera/thermo/IdealGasPhase.h>
#include <cantera/kinetics/InterfaceKinetics.h>
#include <cantera/kinetics/Kinetics.h>
#include <cantera/base/Solution.h>
#include <iostream>
#endif

#ifdef ALAMO_CANTERA
// The actual code is put into a function that can be called from the main program.
void simple_demo()
{
    // Create a new Solution object
  auto sol = Cantera::newSolution("gri30.yaml", "gri30", "None");
    auto gas = sol->thermo();

    // Set the thermodynamic state by specifying T (500 K) P (2 atm) and the mole
    // fractions. Note that the mole fractions do not need to sum to 1.0 - they will
    // be normalized internally. Also, the values for any unspecified species will be
    // set to zero.
    Cantera::AnyMap state;
    state["U"] = 100000.0;
    state["V"] = 1.0e-3;
    state["Y"] = "H2:2.0, O2:2.0";
    
    gas->setState(state);
    std::cout << gas->report() << std::endl;

    Cantera::IdealGasReactor reactor;
    reactor.setSolution(sol);
    Cantera::ReactorNet reactor_net;
    reactor_net.addReactor(reactor);
    double dt = 1e2;
    double t = 0.0;
    double tend = 1e8;
    while (t < tend) {
      t = reactor_net.step();
    }
    //gas->setState_TRX(500.0, 1.0*Cantera::OneAtm, "H2O:1.0, H2:8.0, AR:1.0");

    // Print a summary report of the state of the gas.
    std::cout << gas->report() << std::endl;
}
#endif


// The main program just calls function simple_demo within a 'try' block, and catches
// CanteraError exceptions that might be thrown.
int main()
{
#ifdef ALAMO_CANTERA
    try {
        simple_demo();
    } catch (Cantera::CanteraError& err) {
        std::cout << err.what() << std::endl;
        return 1;
    }
#else
    Util::Abort(INFO,"Needs Cantera");
#endif

    return 0;
}

