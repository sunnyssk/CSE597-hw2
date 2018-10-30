#include "matrix_mpi.h"

template <typename T> int MMat<T>::mpi_rank_;
template <typename T> int MMat<T>::mpi_size_;

template <typename T> MMat<T>::MMat (int rows, int cols) : 
rows_(rows), cols_(cols), slice_rows_((rows - 1) / mpi_size_ + 1), slice_offset_(slice_rows_ * mpi_rank_), buffer_(nullptr) {
    slice_rows_ = ((mpi_rank_ + 1) * slice_rows_ > rows) ? (rows - slice_rows_ * mpi_rank_) : (slice_rows_);
    slice_rows_ = slice_rows_ > 0 ? slice_rows_ : 0;
    buffer_ = new T[slice_rows_ * cols_];
}

template <typename T> T MMat<T>::GetVal (int row, int col) const {
    MPI_Barrier(MPI_COMM_WORLD);
    T retval;
    if(row >= slice_offset_ && row < slice_offset_ + slice_rows_){
        retval = buffer_[];
    }else{
        
    }
    return retval;
}

template <typename T> void MMat<T>::Init (int mpi_size, int mpi_rank) {
    mpi_size_ = mpi_size;
    mpi_rank_ = mpi_rank;
}

template <typename T> void MMat<T>::PrintMatrix (std::ostream output_stream) const {
    MPI_Barrier(MPI_COMM_WORLD);
}