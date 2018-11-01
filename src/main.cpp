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
    
    int rows = 9, cols = 4;
    MMatD A(rows, cols);
    int slice_rows = 0;
    if ((slice_rows = A.SliceRows()) > 0) {
        for(int i = 0; i < slice_rows; i++)
            for(int j = 0; j < cols; j++) A.Elem(i, j) = cols * (A.SliceOffset() + i) + j;
    }
    std::cout << "Rows on process # " << mpi_rank << ": " << A.SliceRows() << std::endl;
    std::cout << "Offset on process # " << mpi_rank << ": " << A.SliceOffset() << std::endl;
    A.Print(std::cout);

    MMatD::Clean();
    MPI_Finalize();
    return 0;
}