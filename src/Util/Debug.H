#ifndef UTIL_DEBUG_H
#define UTIL_DEBUG_H

#include "Util/Util.H"

namespace Util
{
namespace Debug
{

template<class T>
AMREX_FORCE_INLINE
void Probe(
    std::string file, std::string func, int line,
    amrex::FabArray<amrex::BaseFab<T>> &a_mf, int i, int j, int k, int n)
{
    // Loop over boxes owned by this rank
    for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& validbox = mfi.validbox();
        const amrex::Box& box = mfi.grownnodaltilebox();

        if (box.contains(amrex::IntVect(AMREX_D_DECL(i,j,0))))
        {
            const auto& arr = a_mf.const_array(mfi);
            if (validbox.strictly_contains(i,j,k))
                Util::ParallelMessage(file,func,line, "inside: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else if (validbox.contains(i,j,k))
                Util::ParallelMessage(file,func,line, "border: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else
                Util::ParallelMessage(file,func,line, "ghost:  (",i,",",j,",",k,") ",arr(i, j, k, n));
        }
    }
}

AMREX_FORCE_INLINE
int Probe(
    std::string file, std::string func, int line,
    const amrex::MultiFab &a_mf, int _i, int _j, int /*_k*/, int n,
    int pad = 0)
{
    const int t_notexist = 0;
    const int t_inside = 1;
    const int t_border = 2;
    const int t_ghost = 3;
    
    int hit = 0;

    for (int rank = 0; rank < amrex::ParallelDescriptor::NProcs() ; rank++)
    {
        amrex::ParallelDescriptor::Barrier();
        if (rank != amrex::ParallelDescriptor::MyRankInTeam()) continue;

        // Loop over boxes owned by this rank
        for (amrex::MFIter mfi(a_mf,false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& validbox = mfi.validbox();
            const amrex::Box& box = mfi.grownnodaltilebox();
            bool contains = false;

            std::vector<std::vector<Set::Scalar>> vals(2*pad + 1);
            std::vector<std::vector<int>> types(2*pad + 1);
            for (auto & val : vals) val.resize(2*pad + 1,NAN);
            for (auto & type : types) type.resize(2*pad + 1,t_notexist);

            for (int I = 0; I <= 2*pad; I++)
            {
                for (int J = 0; J <= 2*pad; J++)
                {
                    int i = (_i - pad) + I;
                    int j = (_j - pad) + J;
                
                    if (box.contains(amrex::IntVect(AMREX_D_DECL(i,j,0))))
                    {
                        const auto& arr = a_mf.const_array(mfi);
                        vals[I][J] = arr(i,j,0,n);
                        if (validbox.strictly_contains(i,j,0)) types[I][J] = t_inside;
                        else if (validbox.contains(i,j,0))     types[I][J] = t_border;
                        else types[I][J] = t_ghost;
                        contains = true;
                    }
                }
            }

            if (contains)
            {
                hit++;
                std::stringstream ss;
                ss << "  " << std::setw(10) << Color::FG::LightYellow;
                for (int I = 0; I <= 2*pad; I++)
                {
                    int i = (_i - pad) + I;
                    ss << "i=" << std::left << std::setw(13)  << i;
                }
                Util::ParallelMessage(file,func,line,ss.str());
                

                for (int J = 2*pad; J >= 0; J--)
                {
                    int j = (_j - pad) + J;
                
                    std::stringstream ss;
                    ss << Color::FG::LightYellow << "j=" << std::left << std::setw(5) << j << Color::Reset;
                    for (int I = 0; I <= 2*pad; I++)
                    {
                        if (types[I][J] == t_ghost) ss << Color::Dim;
                        if (types[I][J] == t_border) ss << Color::FG::Green;
                        if (types[I][J] == t_inside) ss << Color::FG::Blue;
                        if (types[I][J] == t_notexist) ss << Color::FG::Red;
                        ss << std::left << std::setw(15) << vals[I][J];
                        ss << Color::Reset;
                    }
                    Util::ParallelMessage(file,func,line,ss.str());
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceIntSum(hit);
    return hit;
}


AMREX_FORCE_INLINE
void Probe(
    std::string file, std::string func, int line,
    amrex::iMultiFab &a_mf, int i, int j, int k, int n)
{
    // Loop over boxes owned by this rank
    for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi)
    {
        const amrex::Box& validbox = mfi.validbox();
        const amrex::Box& box = mfi.grownnodaltilebox();

        if (box.contains(amrex::IntVect(AMREX_D_DECL(i,j,0))))
        {
            const auto& arr = a_mf.const_array(mfi);
            if (validbox.strictly_contains(i,j,k))
                Util::ParallelMessage(file,func,line, "inside: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else if (validbox.contains(i,j,k))
                Util::ParallelMessage(file,func,line, "border: (",i,",",j,",",k,") ",arr(i, j, k, n));
            else
                Util::ParallelMessage(file,func,line, "ghost:  (",i,",",j,",",k,") ",arr(i, j, k, n));
        }
    }
}

void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, std::string desc, const amrex::Box domain, int ngrow);

AMREX_FORCE_INLINE
void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, const amrex::Box domain, int ngrow)
{
    Compare(file,func,line,mf,"",domain, ngrow);
}
AMREX_FORCE_INLINE
void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, int ngrow)
{
    Compare(file,func,line,mf,"",amrex::Box::TheUnitBox(), ngrow);
}

AMREX_FORCE_INLINE
void Compare(
    std::string file, std::string func, int line,
    const amrex::MultiFab &mf, std::string desc="")
{    
    Compare(file,func,line,mf,desc,amrex::Box::TheUnitBox(),0);
}

}
}

#endif
