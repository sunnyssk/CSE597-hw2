#ifndef _DEBYE_MPI_H_
#define _DEBYE_MPI_H_

#include <cmath>
#include "mpi_util.h"
#include "matrix_mpi.h"

class MField3D{
public:
    MField3D (int nx, int ny, int nz, double dx, double dy, double dz, int mpi_size, int mpi_rank);
    ~MField3D () { if(buffer_ != nullptr) delete[] buffer_; }

    int Nx () const { return nx_; }
    int Ny () const { return ny_; }
    int Nz () const { return nz_; }
    int N () const { return n_; }
    double Dx () const { return dx_; }
    double Dy () const { return dy_; }
    double Dz () const { return dz_; }
    double * Buffer () const { return buffer_; }

protected:
    int nx_;
    int ny_;
    int nz_;
    int n_;
    int slice_nz_;
    int slice_spacing_;

    double dx_;
    double dy_;
    double dz_;
    double * buffer_;
};

class MDebyeSolver {
public:

    void GenerateSolverMatrix (MField3D const & rhs, double debye_length);
    
protected:
    int mpi_size_;
    int mpi_rank_;

    MMatD* pAmat_;
};

#endif /* _DEBYE_MPI_H_ */