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
    
    int rows = 4, cols = 4;
    MMatD A(rows, cols);
    int slice_rows = 0;
    if((slice_rows = A.SliceRows()) > 0){
        for(int i = 0; i < slice_rows; i++)
            for(int j = 0; j < 4; j++) A(i, j) = 4 * (i + mpi_rank) + j;
    }
    if(mpi_rank == 0){
        double * res = new double[rows * cols];
        for(int i = 1; i < mpi_size; i++) MPI_Recv(A.Buffer() + i * A.SliceRows() * A.Cols(), (rows % A.SliceRows()) * A.Cols(), MPI_DOUBLE, i, i, MPI_COMM_WORLD, nullptr);
        std::cout << "The matrix is:\n";
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++) std::cout << res[i * cols + j] << " ";
            std::cout << "\n";
        }
        delete[] res;
    }else{
        MPI_Send(A.ReadOnlyBuffer(), A.SliceRows() * A.Cols(), MPI_DOUBLE, 0, mpi_rank, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}