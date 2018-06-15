#include "WriteMetaData.H"

void IO::WriteMetaData(std::string plot_file) 
{
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

  auto starttime = std::chrono::system_clock::now();
  std::time_t now = std::chrono::system_clock::to_time_t(starttime);
  metadatafile << "# RUN DETAILS" << std::endl;
  metadatafile << "# ===========" << std::endl;
  metadatafile << "Simulation_start_time = " << std::ctime(&now);
  metadatafile << "Number_of_processors = " << amrex::ParallelDescriptor::NProcs() << std::endl;
  metadatafile << std::endl;

  metadatafile << "# PARAMETERS" << std::endl;
  metadatafile << "# ==========" << std::endl;
  amrex::ParmParse pp;
  pp.dumpTable(metadatafile,true);
  
  metadatafile.close();
}

