#include "WriteMetaData.H"
#include <sstream>
#include <functional>
#include "Util/Util.H"

/// \todo We need to update ParmParse so that the FORMATTED input file gets recorded permanently.

namespace IO
{

// GLOBAL VARIABLES
unsigned long hash = 0;
std::chrono::time_point<std::chrono::system_clock> starttime_cr;
std::time_t starttime = 0;
int percent = -1;

void WriteMetaData(std::string plot_file, Status status, int per) 
{
    if (amrex::ParallelDescriptor::IOProcessor())
        {
            amrex::ParmParse pp;

            std::chrono::time_point<std::chrono::system_clock> now_cr = std::chrono::system_clock::now();
            std::time_t now = std::chrono::system_clock::to_time_t(now_cr);

            if (!starttime)
            {
                starttime_cr = now_cr;
                starttime = now;
            } 

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
            {
                char buffer [80];
                std::strftime(buffer,80,"%Y-%m-%d %H:%M:%S",std::localtime(&starttime));
                std::string timefmt(buffer);
                metadatafile << "Simulation_start_time = " << timefmt << std::endl;
            }
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
                    char buffer[80];
                    std::strftime(buffer,80,"%Y-%m-%d %H:%M:%S",std::localtime(&now));
                    std::string timefmt(buffer);
                    metadatafile << "Simulation_end_time = " << timefmt << std::endl;
                }

            auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(now_cr - starttime_cr);
            metadatafile << "Simulation_run_time = " << (float)milliseconds.count()/1000.0 << " " << std::endl;

            #ifdef GIT_DIFF_OUTPUT
            {
                std::ifstream src(GIT_DIFF_OUTPUT,std::ios::binary);
                std::ofstream dst(plot_file+"/diff.html",std::ios::binary);
                dst << src.rdbuf();
                src.close();
                dst.close();
            }
            #endif
            {
                std::ifstream src(GIT_DIFF_PATCH_OUTPUT,std::ios::binary);
                if (src.is_open())
                {
                    std::ofstream dst(plot_file+"/diff.patch",std::ios::binary);
                    dst << src.rdbuf();
                    src.close();
                    dst.close();
                }
                else
                {
                    Util::Warning(INFO,"Could not open ",GIT_DIFF_PATCH_OUTPUT);
                }
                
            }

            metadatafile.close();
        }
}

}
