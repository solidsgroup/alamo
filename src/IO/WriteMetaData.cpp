#include "WriteMetaData.H"
#include <sstream>
#include <functional>

/// \todo We need to update ParmParse so that the FORMATTED input file gets recorded permanently.

namespace IO
{

// GLOBAL VARIABLES
unsigned long hash = 0;
std::time_t starttime = 0;
int percent = -1;

void WriteMetaData(std::string plot_file, Status status, int per) 
{
	if (amrex::ParallelDescriptor::IOProcessor())
		{
			amrex::ParmParse pp;

			std::time_t now = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());

			if (!starttime) starttime = now;

			if (!hash)
				{
					std::stringstream ss;
					pp.dumpTable(ss);
					std::string dumptable = ss.str();
					hash =  std::hash<std::string>{}(dumptable + std::ctime(&starttime));
				}

			if (per > percent) percent = per;

			std::ofstream metadatafile;
			metadatafile.open(plot_file+"/metadata",std::ios_base::out);
			metadatafile << "# COMPILATION DETAILS" << std::endl;
			metadatafile << "# ===================" << std::endl;
			metadatafile << "Git_commit_hash = " << METADATA_GITHASH << std::endl;
			metadatafile << "AMReX_version = " << amrex::Version() << std::endl;
			metadatafile << "Dimension = " << AMREX_SPACEDIM << std::endl;
			metadatafile << "User = " << METADATA_USER  << std::endl;
			metadatafile << "Platform = " << METADATA_PLATFORM  << std::endl;
			metadatafile << "Compiler = " << METADATA_COMPILER  << std::endl;
			metadatafile << "Compilation_Date = " << METADATA_DATE  << std::endl;
			metadatafile << "Compilation_Time = " << METADATA_TIME  << std::endl;
			metadatafile << std::endl;

			metadatafile << "# PARAMETERS" << std::endl;
			metadatafile << "# ==========" << std::endl;

			pp.dumpTable(metadatafile,true);
  
			metadatafile << "# RUN DETAILS" << std::endl;
			metadatafile << "# ===========" << std::endl;
			metadatafile << "HASH = " << hash << std::endl;
			metadatafile << "Simulation_start_time = " << std::ctime(&starttime);
			metadatafile << "Number_of_processors = " << amrex::ParallelDescriptor::NProcs() << std::endl;
			if (status == Status::Running)        metadatafile << "Status = Running";
			else if (status == Status::Complete)  metadatafile << "Status = Complete";
			else if (status == Status::Abort)     metadatafile << "Status = Abort";
			else if (status == Status::Segfault)  metadatafile << "Status = Segfault";
			else if (status == Status::Interrupt) metadatafile << "Status = Interrupt";

			if (percent>0) metadatafile << " (" << percent << "%)";
			metadatafile << std::endl;

			if (status != Status::Running)
				{
					metadatafile << "Simulation_end_time = " << std::ctime(&starttime);
				}
			metadatafile << "Simulation_run_time = " << std::difftime(now,starttime) << " seconds " << std::endl;

			#ifdef GIT_DIFF_OUTPUT
			std::ifstream src(GIT_DIFF_OUTPUT,std::ios::binary);
			std::ofstream dst(plot_file+"/diff.html",std::ios::binary);
			dst << src.rdbuf();
			src.close();
			dst.close();
			#endif

			metadatafile.close();
		}
}

}
