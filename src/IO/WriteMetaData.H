//
// This provide an IO routine for generating run-specific output.
// Every Alamo run produces a metadata file and (if applicable) a diff.patch
// file that reflects the exact state when the simulation was run.
// Note that the file :code:`WriteMetaData.cpp` file is _always_ recompiled
// upon :code:`make`, in order to pull the latest information about the local
// state at build.
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
