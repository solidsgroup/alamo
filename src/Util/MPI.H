#ifndef UTIL_MPI_H
#define UTIL_MPI_H

#include <vector>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_ParallelContext.H>


namespace Util
{
namespace MPI
{



template <class T>
int Allgather(std::vector<T>& a_data)
{
    // Gather information about how many sites were found on each processor
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int my_num = a_data.size();
    int num = my_num;

    // Communicate the total array size to everyone
    amrex::ParallelAllReduce::Sum(num, MPI_COMM_WORLD);
    // Temporary buffers to receive data from all procs
    std::vector<T> a_data_all(num);
    // Send information about how many sites on each proc to all procs
    std::vector<int> nsites_procs(nprocs);
    MPI_Allgather(  &my_num, 1, amrex::ParallelDescriptor::Mpi_typemap<int>::type(),
                    nsites_procs.data(), 1, amrex::ParallelDescriptor::Mpi_typemap<int>::type(),
                    MPI_COMM_WORLD);
    // Calculate the offset for each
    std::vector<int> nsites_disp(nprocs);
    for (int i = 0; i < nprocs; i++)
    {
        nsites_disp[i] = 0;
        for (int j = 0; j < i; j++) nsites_disp[i] += nsites_procs[j];
    }
    // Store the MPI datatype for each
    MPI_Datatype mpi_type = amrex::ParallelDescriptor::Mpi_typemap<T>::type();
    MPI_Allgatherv(
        a_data.data(), my_num, mpi_type,
        a_data_all.data(), nsites_procs.data(), nsites_disp.data(), mpi_type,
        MPI_COMM_WORLD);
    // Swap out the data so the buffers are no longer needed.
    a_data.swap(a_data_all);
    a_data_all.clear();
    return 0;
}

inline void Bcast(  std::string &s,
                    int root = amrex::ParallelDescriptor::IOProcessorNumber())
{
    using namespace amrex;

    int me = ParallelDescriptor::MyProc();

    // --- Broadcast length ---
    int len = (me == root) ? static_cast<int>(s.size()) : 0;
    ParallelDescriptor::Bcast(&len, 1, root);

    // --- Resize on receivers ---
    if (me != root)
        s.resize(len);

    // --- Broadcast character data ---
    if (len > 0)
        ParallelDescriptor::Bcast(s.data(), len, root);
}





inline void AllGather(  const std::string &local,
                        std::vector<std::string> &out)
{
    using namespace amrex;

    int nprocs = ParallelDescriptor::NProcs();
    int me     = ParallelDescriptor::MyProc();

    // Step 1: gather lengths from all ranks
    int local_len = static_cast<int>(local.size());
    std::vector<int> lengths(nprocs);

    for (int p = 0; p < nprocs; ++p) {
        int len = (me == p) ? local_len : 0;
        ParallelDescriptor::Bcast(&len, 1, p);
        lengths[p] = len;
    }

    // Step 2: compute offsets
    std::vector<int> offsets(nprocs, 0);
    int total_len = lengths[0];
    for (int i = 1; i < nprocs; ++i) {
        offsets[i] = offsets[i-1] + lengths[i-1];
        total_len += lengths[i];
    }

    // Step 3: prepare local buffer
    std::vector<char> local_buf(local_len);
    if (local_len > 0)
        std::memcpy(local_buf.data(), local.data(), local_len);

    // Step 4: broadcast character buffers from each rank
    std::vector<char> all_buf(total_len);
    for (int p = 0; p < nprocs; ++p) {
        int len = lengths[p];
        if (me == p) {
            // my buffer is already local_buf
            if (len > 0)
                std::memcpy(all_buf.data() + offsets[p], local_buf.data(), len);
        }
        // broadcast this rank's characters to everyone
        ParallelDescriptor::Bcast(all_buf.data() + offsets[p], len, p);
    }

    // Step 5: unpack strings
    out.resize(nprocs);
    for (int i = 0; i < nprocs; ++i) {
        out[i] = std::string(all_buf.data() + offsets[i], lengths[i]);
    }
}




}
}




#endif
