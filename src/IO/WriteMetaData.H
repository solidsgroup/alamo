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
