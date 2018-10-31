#include "debye_mpi.h"

MField3D (int nx, int ny, int nz, double dx, double dy, double dz, int mpi_size, int mpi_rank) :
nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz), mpi_size_(mpi_size), mpi_rank_(mpi_rank) {

}

void MDebyeSolver::GenerateSolverMatrix (MField3D const & rhs, double debye_length) {
    int nx = rhs.Nx(), ny = rhs.Ny(), nz = rhs.Nz();
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nltrim * nztrim;
    double dx = rhs.Dx(), dy = rhs.Dy(), dz = rhs.Dz();
    double dx2 = dx * dx, dy2 = dy * dy, dz2 = dz * dz, ldi2 = 1 / (debye_length * debye_length);
    if (pAmat_) delete pAmat_;
    pAmat_ = new MMatD(ntrim, ntrim);
    MMatD &A = *pAmat_;

    // TODO: modify k; avoid modifying other elements
    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++) {
        for (int j = 1; j <= nytrim; j++) {
            for (int k = 1; k <= nztrim; k++) {
                A(TRI(i, j, k), TRI(i, j, k)) = -(2 / dx2 + 2 / dy2 + 2 / dz2 + ldi2);
                if(i > 1) A(TRI(i, j, k), TRI(i - 1, j, k)) = 1 / dx2;
                if(j > 1) A(TRI(i, j, k), TRI(i, j - 1, k)) = 1 / dy2;
                if(k > 1) A(TRI(i, j, k), TRI(i, j, k - 1)) = 1 / dz2;
                if(i < nxtrim) A(TRI(i, j, k), TRI(i + 1, j, k)) = 1 / dx2;
                if(j < nytrim) A(TRI(i, j, k), TRI(i, j + 1, k)) = 1 / dy2;
                if(k < nztrim) A(TRI(i, j, k), TRI(i, j, k + 1)) = 1 / dz2;
            }
        }
    }
    #undef TRI
}