#ifndef _MATRIX_MPI_H_
#define _MATRIX_MPI_H_

#include "mpi_util.h"

// Matrices are stored piecewise in processes, with (slice_size_) rows each piece.
// Data is stored with column number changing in the most rapid manner.

template <typename T> class MMat {
public:
    MMat (int rows, int cols) : rows_(rows),
                                cols_(cols),
                                slice_rows_(rows / mpi_size_),
                                slice_offset_(slice_rows_ * mpi_rank_),
                                buffer_(nullptr) {
                                    slice_rows_ = rows % slice_rows_;
                                    buffer_ = new T[slice_rows_ * cols_];
                                }

    ~MMat () { delete[] buffer_; }

    static void Init (int mpi_size, int mpi_rank);

    int Rows() const { return rows_; }
    int Cols() const { return cols_; }
    int SliceRows() const { return slice_rows_; }
    int SliceOffset() const { return slice_offset_; }
    T const * const ReadOnlyBuffer() const { return buffer_; }
    T * const Buffer() { return buffer_; }

    T& operator () (int slice_row, int col) { return buffer_[slice_row * cols_ + col]; }
    const T& operator () (int slice_row, int col) const { return buffer_[slice_row * cols_ + col]; }


protected:
    static int mpi_size_;
    static int mpi_rank_;
    int rows_;
    int cols_;
    int slice_rows_;
    int slice_offset_;
    T* buffer_;
};

template class MMat<int>;
template class MMat<float>;
template class MMat<double>;

typedef MMat<int> MMatI;
typedef MMat<float> MMatF;
typedef MMat<double> MMatD;

#endif /* _MATRIX_MPI_H_ */