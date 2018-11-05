#include "mpi_util.h"

void MPIArraySync (double * array, int size, int mpi_size, int mpi_rank) {
    int mpi_interval = (size - 1) / mpi_size + 1, slice_rows = 0, offset = 0;
    for (int i = 0; i < mpi_size; i++) {
        slice_rows = ((i + 1) * mpi_interval > size) ? (size - mpi_interval * i) : (mpi_interval);
        slice_rows = slice_rows > 0 ? slice_rows : 0;
        if(slice_rows > 0){
            MPI_Bcast(array + offset, slice_rows, MPI_DOUBLE, i, MPI_COMM_WORLD);
        }
        offset += mpi_interval;
    }
}

void MPIArraySyncAllGatherv (double * src, double * tar, int size, int mpi_size, int mpi_rank) {
    int mpi_interval = (size - 1) / mpi_size + 1, slice_rows = 0;
    int *count = new int[mpi_size], *offset = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        slice_rows = ((i + 1) * mpi_interval > size) ? (size - mpi_interval * i) : (mpi_interval);
        slice_rows = slice_rows > 0 ? slice_rows : 0;
        count[i] = slice_rows;
        offset[i] = (i > 0) ? (offset[i - 1] + count[i - 1]) : 0;
    }
    MPI_Allgatherv(src + offset[mpi_rank], count[mpi_rank], MPI_DOUBLE, tar, count, offset, MPI_DOUBLE, MPI_COMM_WORLD);
    delete[] count;
    delete[] offset;
}

void MPIGatherMax (double * target, int mpi_size, int mpi_rank) {
    double buffer = 0.0;
    for (int i = 0; i < mpi_size; i++) {
        if(i == mpi_rank) buffer = *target;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&buffer, 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
        if(buffer > *target) *target = buffer;
    }
}

void MPIStart (int & mpi_size, int & mpi_rank) {
    MPI_Init(nullptr, nullptr);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
}
