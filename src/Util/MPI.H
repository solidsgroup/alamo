#ifndef UTIL_MPI_H
#define UTIL_MPI_H

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
    MPI_Allgather(&my_num, 1, amrex::ParallelDescriptor::Mpi_typemap<int>::type(),
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


}
}



#endif
