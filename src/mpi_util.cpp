#include "mpi_util.h"

void MPIStart (int & mpi_size, int & mpi_rank) {
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
}
