#ifndef _MATRIX_MPI_H_
#define _MATRIX_MPI_H_

#include "mpi_util.h"

template <typename T> class MMat {
public:
    MMat (int rows, int cols) : 

    static void Init (int mpi_size, int mpi_rank);

    int Rows() const { return rows_; }
    int Cols() const { return cols_; }

protected:
    static int mpi_size_;
    static int mpi_rank_;
    int rows_;
    int cols_;
    int slice_size_;
    T* buffer_;
};

template class MMat<int>;
template class MMat<float>;
template class MMat<double>;

typedef MMat<int> MMatI;
typedef MMat<int> MMatF;
typedef MMat<int> MMatD;

#endif /* _MATRIX_MPI_H_ */