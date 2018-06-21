#include <string>

#include "AMReX.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

namespace Util
{

std::string GetFileName();

void SignalHandler(int s);

void Initialize (int argc, char* argv[]);

void Finalize ();

void Terminate(const char * msg, int signal, bool backtrace);

void Abort (const std::string& msg);
void Abort (const char * msg);

}
