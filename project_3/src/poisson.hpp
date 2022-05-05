#include "mpi.h"
#include "DistributedMatrix.hpp"

#ifndef poissonHEADERDEF
#define poissonHEADERDEF

void jacobi_iteration(const DistributedMatrix &a, const DistributedMatrix &f, DistributedMatrix &b);
double sum_squares(const DistributedMatrix &a, const DistributedMatrix &b);

#endif
