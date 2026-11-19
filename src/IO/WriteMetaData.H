//
// This provide an IO routine for generating run-specific output.
// Every Alamo run records build metadata and copies the generated Git patches.
// The build refreshes the patches separately and recompiles the metadata object
// when the development merge base changes so the patch and recorded hash agree.
//
// :bdg-warning-line:`This is standard infrastructure code; make sure you know what you ard doing before you change it.`
//

#ifndef IO_WRITEMETADATA_H
#define IO_WRITEMETADATA_H

#include <ostream>
#include <iostream>
#include <string>
#include <ctime>
#include <fstream>
#include <chrono>

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

namespace IO
{

enum Status{
    Running,
    Complete,
    Segfault,
    Interrupt,
    Abort
};

    void WriteMetaData(std::string plot_file, Status status = Status::Running, int percent = -1);


}

#endif
