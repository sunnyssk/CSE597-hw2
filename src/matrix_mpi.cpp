#include "matrix_mpi.h"

template <typename T> int MMat<T>::mpi_rank_;
template <typename T> int MMat<T>::mpi_size_;

template <typename T> void MMat<T>::Init (int mpi_size, int mpi_rank) {
    mpi_size_ = mpi_size;
    mpi_rank_ = mpi_rank;
}

// template <typename T> T& Val (int clice_row, int col) {
//     if(row >= slice_offset_ && row < slice_offset_ + slice_size_){
//
//     }else{
//
//     }
// }