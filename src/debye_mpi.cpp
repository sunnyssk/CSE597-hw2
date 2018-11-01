#include "debye_mpi.h"

MField3D::MField3D (int nx, int ny, int nz, double dx, double dy, double dz, int mpi_size, int mpi_rank) :
nx_(nx), ny_(ny), nz_(nz), dx_(dx), dy_(dy), dz_(dz) {

}

void MDebyeSolver::GenerateSolverMatrix (MField3D const & rhs, double debye_length) {
    int nx = rhs.Nx(), ny = rhs.Ny(), nz = rhs.Nz();
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nltrim * nztrim;
    double dx = rhs.Dx(), dy = rhs.Dy(), dz = rhs.Dz();
    double dx2 = dx * dx, dy2 = dy * dy, dz2 = dz * dz, ldi2 = 1 / (debye_length * debye_length);
    if (pAmat_) delete pAmat_;
    pAmat_ = new MMatD(ntrim, ntrim);
    MMatD &A = *pAmat_;
    int row_offset = A.SliceOffset(), row_terminate = row_offset + A.SliceRows();

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int k = 1; k <= nztrim; k++) {
        if (TRI(0, 0, k + 1) <= row_offset || TRI(0, 0, k) > row_terminate) continue;
        for (int j = 1; j <= nytrim; j++) {
            if (TRI(0, j + 1, k) <= row_offset || TRI(0, j, k) > row_terminate) continue;
            for (int i = 1; i <= nxtrim; i++) {
                if (TRI(i + 1, j, k) <= row_offset || TRI(i, j, k) > row_terminate) continue;
                A.Elem(TRI(i, j, k) - row_offset, TRI(i, j, k)) = -(2 / dx2 + 2 / dy2 + 2 / dz2 + ldi2);
                if(i > 1) A.Elem(TRI(i, j, k) - row_offset, TRI(i - 1, j, k)) = 1 / dx2;
                if(j > 1) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j - 1, k)) = 1 / dy2;
                if(k > 1) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j, k - 1)) = 1 / dz2;
                if(i < nxtrim) A.Elem(TRI(i, j, k) - row_offset, TRI(i + 1, j, k)) = 1 / dx2;
                if(j < nytrim) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j + 1, k)) = 1 / dy2;
                if(k < nztrim) A.Elem(TRI(i, j, k) - row_offset, TRI(i, j, k + 1)) = 1 / dz2;
            }
        }
    }
    #undef TRI
}

int MDebyeSolver::JacobiIterativeSolve (double err_threshold, MField3D & res_container, double * iter_err_array) {
    int nx = res_container.Nx(), ny = res_container.Ny(), nz = res_container.Nz();
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    int iter_cnt = 0;
    int slice_rows = pAmat_->SliceRows(), slice_offset = pAmat->SliceOffset();

    double *aux_vector1 = new double[ntrim], *aux_vector2 = new double[ntrim];
    double *x_prev = aux_vector1, *x_next = aux_vector2;
    MatD *b = pfield_;
    double errmax = 0.0;

    // Assign the initial guess to the input vector
    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) (*x_prev)(TRI(i, j, k)) = res_container.GetVal(i, j, k);
    #undef TRI

    // Iterative solver
    do {
        errmax = 0.0;
        if (iter_cnt++ >= MAX_ITER_NUM - 1) {
            std::cout << "Iteration number exceeds limit! Exit automatically." << std::endl;
            break;
        }
        for (int i = slice_offset; i < slice_offset + slice_rows; i++) {
            x_next[i] = pfield_[i];
            for (int j = 0; j < i; j++) x_next[i] -= pAmat_->Elem(i, j) * x_prev[j];
            for (int j = i + 1; j < ntrim; j++) x_next[i] -= pAmat_->Elem(i, j) * x_prev[j];
            x_next[i] /= pAmat_->Elem(i, i);
            double newerr = fabs(x_next[i] - x_prev[i]);
            errmax = errmax > newerr ? errmax : newerr;
        }
        // TODO: COMMUNICATION ON RESULT VECTOR & ERRMAX
        for (int i = 0; i < mpi_rank_; i++)
        double *ptmp = x_prev;                                                  // swap two vectors
        x_prev = x_next;
        x_next = ptmp;
        std::cout << "Iteration round #" << iter_cnt << ", maximum error: " << errmax << "\n";
        iter_err_array[iter_cnt - 1] = errmax;
    } while (errmax >= err_threshold || &aux_vector1 != x_next);                // only exits iteration when the original field is updated

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) res_container(i, j, k) = (*x_next)(TRI(i, j, k));
    #undef TRI

    char *buffer = new char[256];
    MemSizeOutput(buffer);
    delete[] buffer;
    delete[] aux_vector1;
    delete[] aux_vector2;
    return iter_cnt;
}

void MDebyeSolver::RhsInput (MField3D const & field) {
    int nx = field.Nx(), ny = field.Ny(), nz = field.Nz();
    int nxtrim = nx - 2, nytrim = ny - 2, nztrim = nz - 2, nltrim = nxtrim * nytrim, ntrim = nxtrim * nytrim * nztrim;
    if (pfield_) delete pfield_;
    pfield_ = new double(ntrim, 1);

    #define TRI(i, j, k) (((k) - 1) * nltrim + ((j) - 1) * nxtrim + (i) - 1)    // trimmed indices
    for (int i = 1; i <= nxtrim; i++)
        for (int j = 1; j <= nytrim; j++)
            for (int k = 1; k <= nztrim; k++) pfield_[TRI(i, j, k)] = field(i, j, k);
    #undef TRI
}