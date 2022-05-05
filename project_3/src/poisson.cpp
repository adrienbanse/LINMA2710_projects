#include "poisson.hpp"

void jacobi_iteration(const DistributedMatrix &a, const DistributedMatrix &f, DistributedMatrix &b){
    int n, i, j;
    double h;
    n = a.Size(1) - 1;
    h = 1 / n;
    for (i = a.FirstWriteRow(); i <= a.LastWriteRow(); i++){
        for (j = a.FirstWriteColumn(); j <= a.LastWriteColumn(); j++){
            if (i == 1 || i == n + 1 || j == 1 || j == n + 1)
                b(i, j) = a.Read(i, j);
            else
                b(i, j) = (a.Read(i - 1, j) + a.Read(i + 1, j) + a.Read(i, j - 1) + a.Read(i, j + 1) - h * h * f.Read(i, j)) / 4;
        }
    }
}

double sum_squares(const DistributedMatrix &a, const DistributedMatrix &b){
    double delta, block_sum, sum;
    int i, j;
    block_sum = 0.;
    for (i = a.FirstWriteRow(); i <= a.LastWriteRow(); i++){
        for (j = a.FirstWriteColumn(); j <= a.LastWriteColumn(); j++){
            delta = b.Read(i, j) - a.Read(i, j);
            block_sum += delta * delta;
        }
    }
    MPI_Allreduce(&block_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sum;
}