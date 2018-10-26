#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <mpi.h>
#include "mpi_util.h"
#include "matrix_mpi.h"

int main (int argc, char** argv) {
    int mpi_size = 0, mpi_rank = 0;
    MPIStart(mpi_size, mpi_rank);
    MMatD::Init(mpi_size, mpi_rank);
    
    MMatD A(3, 3);

    MPI_Finalize();
    return 0;
}