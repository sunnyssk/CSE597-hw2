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
    if((slice_rows = A.SliceRows()) > 0){
        for(int i = 0; i < slice_rows; i++)
            for(int j = 0; j < cols; j++) A(i, j) = cols * (A.SliceOffset() + i) + j;
    }
    std::cout << "Rows on process # " << mpi_rank << ": " << A.SliceRows() << std::endl;
    if(mpi_rank == 0){
        double * res = new double[rows * cols];
        memcpy(res, A.Buffer(), sizeof(double) * A.SliceRows() * A.Cols());
        int restrows = rows;
        for(int i = 1; i < mpi_size; i++){
            restrows -= A.SliceRows();
            int srcrows = (restrows >= A.SliceRows()) ? A.SliceRows() : restrows;
            if(srcrows > 0) MPI_Recv(res + i * A.SliceRows() * A.Cols(), srcrows * A.Cols(), MPI_DOUBLE, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        std::cout << "The matrix is:\n";
        for(int i = 0; i < rows; i++){
            for(int j = 0; j < cols; j++) std::cout << res[i * cols + j] << " ";
            std::cout << "\n";
        }
        delete[] res;
    }else{
        if(A.SliceRows() > 0) MPI_Send(A.Buffer(), A.SliceRows() * A.Cols(), MPI_DOUBLE, 0, mpi_rank, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}