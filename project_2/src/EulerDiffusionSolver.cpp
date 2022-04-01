#include <cmath>
#include <cstdio>
#include <cassert>
#include <vector>
#include <iostream>
#include "globals.hpp"
#include "EulerDiffusionSolver.hpp"

#define I(a, b) (b)*(nx+1)+(a)

SparseMatrix* EulerDiffusionSolver::computeAb(Vector &b){
    int M = 5;
    int d = (nx + 1) * (ny + 1);
    int *rowidx = new int[d * M];
    double *nzval = new double[d * M];
    int *counterCol = (int *) calloc(d, sizeof(int));
    int nnz = 0;

    int i, j, col;
    double x, y, D;
    for (i = 0; i <= nx; i++){
        for (j = 0; j <= ny; j++){
            x = i * dx;
            y = j * dy;

            if (i == 0) b[I(i, j)] = p.getBoundaryx0(y);
            else if (i == nx) b[I(i, j)] = p.getBoundaryx1(y);
            else if (j == 0) b[I(i, j)] = p.getBoundaryy0(x);
            else if (j == ny) b[I(i, j)] = p.getBoundaryy1(x);
            else{
                b[I(i, j)] = dt * p.getf(x, y);
                D = p.getD(x, y);

                // add (alphax * D(x_i, y_j)) at (I(i, j), I(i-1, j)) and (I(i, j), I(i+1, j))
                col = I(i - 1, j);
                rowidx[col * M + counterCol[col]] = I(i, j);
                nzval[col * M + counterCol[col]] = alphax * D;
                counterCol[col]++;

                col = I(i + 1, j);
                rowidx[col * M + counterCol[col]] = I(i, j);
                nzval[col * M + counterCol[col]] = alphax * D;
                counterCol[col]++;

                // add (alphay * D(x_i, y_j)) at (I(i, j), I(i, j-1)) and (I(i, j), I(i, j+1))
                col = I(i, j - 1);
                rowidx[col * M + counterCol[col]] = I(i, j);
                nzval[col * M + counterCol[col]] = alphay * D;
                counterCol[col]++;

                col = I(i, j + 1);
                rowidx[col * M + counterCol[col]] = I(i, j);
                nzval[col * M + counterCol[col]] = alphay * D;
                counterCol[col]++;

                // add (1 - 2 D(x_i, y_j)[alphax + alphay]) at (I(i, j), I(i, j))
                col = I(i, j);
                rowidx[col * M + counterCol[col]] = I(i, j);
                nzval[col * M + counterCol[col]] = 1 - 2 * D * (alphax + alphay);
                counterCol[col]++;

                nnz += 5;
            }   
        }
    }

    SparseMatrix *Ap = new SparseMatrix(d, d, nnz, M, rowidx, nzval);

    delete[] rowidx;
    delete[] nzval;
    free(counterCol);

    return Ap;
}

EulerDiffusionSolver::EulerDiffusionSolver(int nt, int nx, int ny, const DiffusionProblem &p):p(p){
    this->nt = nt;
    this->nx = nx;
    this->ny = ny;

    dt = p.getLt() / nt;
    dx = p.getLx() / nx;
    dy = p.getLy() / ny;

    alphax = dt / (dx * dx);
    alphay = dt / (dy * dy);

    b = new Vector((nx + 1) * (ny + 1));
    A = computeAb(*(b));
}

void EulerDiffusionSolver::solve(Vector **sol){
    int i, j, k;
    double x, y;
    Vector sol0((nx + 1) * (ny + 1));

    // intial condition
    for (i = 0; i <= nx; i++){
        for (j = 0; j <= ny; j++){
            x = i * dx;
            y = j * dy;
            sol0[I(i, j)] = p.getInitial(x, y);
        }
    }
    *(sol[0]) = sol0;

    // temporal iterations
    for (k = 1; k <= nt; k++){
        *(sol[k]) = *(A) * *(sol[k - 1]) + *(b);
        printf("%d / %d\n", k, nt);
    }

    // free A and b
    delete A;
    delete b;
}
