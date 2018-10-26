#include "matrix_mpi.h"

template <typename T> void MMat<T>::Init (int mpi_size, int mpi_rank) {
    mpi_size_ = mpi_size;
    mpi_rank_ = mpi_rank;
}