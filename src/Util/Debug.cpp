#include "Util/Debug.H"
#include "IO/ParmParse.H"
#include "Util/MPI.H"
#include "IO/FileNameParse.H"
#include "Util/Util.H"

#include <filesystem>


// This code is not included in coverage count since it is for debugging purposes only
// LCOV_EXCL_START

namespace Util
{
namespace Debug
{

bool Compare_on = false;
std::map<std::string,int>  Compare_map;
int Compare_cnt = 0;

void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &a_mf, std::string desc,
    amrex::Box domain, int /*ngrow*/)
{
    IO::ParmParse pp;
    // Switch to enableor disable debug comparisons.
    pp.query_default("util.debug.compare",Compare_on,false);
    if (!Compare_on) return;

    Compare_cnt++;

    Set::Scalar tolerance = 1E-10;

    std::ostringstream messagestream;
    std::string filesan = file;
    Util::String::ReplaceAll(filesan,"/","_");
    messagestream << filesan << "_" << func << "_" << line << desc;
    std::string key = messagestream.str();

    int cntr = 0;
    if (Compare_map.find(key) == Compare_map.end())
        Compare_map[key] = 0;
    else
        Compare_map[key] ++;
    cntr = Compare_map[key];

    messagestream << "_" << cntr;

    std::string name = messagestream.str();
    std::vector<std::string> allnames;
    Util::MPI::Allgather(name,allnames);
    for (auto & othername : allnames)
    {
        if (name != othername)
            Util::ParallelAbort(file,func,line, "MPI Paths have diverged: I am trying to read ", 
                                name, " but some other process is trying to read ", 
                                othername);
    }

    amrex::BoxArray ba = a_mf.boxArray();
    amrex::DistributionMapping dm = a_mf.DistributionMap();
    int nghost = 2;
    int ncomp = a_mf.nComp();
    if (domain != amrex::Box::TheUnitBox())
    {
        for (auto &b : ba.boxList()) b = b.grow(2) & domain;
        nghost = 0;
    }
    amrex::MultiFab mf(ba,dm,ncomp,nghost);
    amrex::MultiFab::Copy(mf,a_mf,0,0,ncomp,0);
    if (amrex::ParallelDescriptor::NProcs() == 1)
    {
        Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ","Storing ",name);
        Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ",name, " - ncomp ",ncomp);
        std::filesystem::create_directory("checks");
        amrex::VisMF::Write(mf, "checks/" + name);
    }
    else
    {
        Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ","Comparing ",name);
        {
            amrex::MultiFab mftmp;
            amrex::VisMF::Read(mftmp,"checks/" + name);
            if (ba != mftmp.boxArray())
            {
                Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ","our boxarray ",ba);
                Util::ParallelMessage(file,func,line,"[",Compare_cnt,"] ","saved boxarray ",mftmp.boxArray());
                Util::Warning(file,func,line,"[",Compare_cnt,"] ","different box arrays !! !! !!");
            }
        }

        amrex::MultiFab mforig(ba,mf.distributionMap,ncomp,nghost);
        amrex::VisMF::Read(mforig,"checks/" + name);

        amrex::MultiFab mfdiff(ba, mf.distributionMap, ncomp, nghost);
        amrex::MultiFab::Copy(mfdiff, mf,0,0,ncomp,nghost);
        amrex::MultiFab::Subtract(mfdiff,mforig,0,0,ncomp,nghost);
        mfdiff.abs(0,ncomp,nghost);

        for (int n = 0 ; n < ncomp; n++)
        {
            Set::Scalar max = mfdiff.max(n,nghost);
            auto argmax = mfdiff.maxIndex(n,nghost);
            if (max > tolerance)
            {
                amrex::ParallelDescriptor::Barrier();
                Util::ParallelMessage(file,func,line,"ours (n=",n,")");
                amrex::ParallelDescriptor::Barrier();
                Util::Debug::Probe(file,func,line,mf,argmax[0],argmax[1],0,n,3);

                amrex::ParallelDescriptor::Barrier();
                Util::ParallelMessage(file,func,line,"original (n=",n,")");
                amrex::ParallelDescriptor::Barrier();
                Util::Debug::Probe(file,func,line,mforig,argmax[0],argmax[1],0,n,3);

                amrex::ParallelDescriptor::Barrier();
                Util::ParallelMessage(file,func,line,"diff (n=",n,")");
                amrex::ParallelDescriptor::Barrier();
                Util::Debug::Probe(file,func,line,mfdiff,argmax[0],argmax[1],0,n,3);

                messagestream << "_diff";
                amrex::VisMF::Write(mfdiff,messagestream.str());
                Util::Warning(file,func,line,"[",Compare_cnt,"] ",name);
                Util::ParallelAbort(file,func,line,"[",Compare_cnt,"] ",max," at ", argmax);
            }
        }
    }
    //if (Compare_cnt == 169) Util::ParallelAbort(INFO,"exiting");
}
}
}
// LCOV_EXCL_START
