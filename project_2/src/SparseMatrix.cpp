#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(int numRows, int numCols)
{
    assert(numRows > 0);
    assert(numCols > 0);
    m = numRows;
    n = numCols;
    M = 0;
    nnz = 0;
    rowidx = new int[0];
    nzval = new double[0];
}

SparseMatrix::SparseMatrix(int m, int n, int nnz, int M, int* rowidx, double* nzval)
{
    assert(m > 0);
    assert(n > 0);
    this->m = m;
    this->n = n;
    this->M = M;
    this->nnz = nnz;
    this->rowidx = new int[n * M];
    this->nzval = new double[n * M];
    memcpy(this->rowidx, rowidx, n * M * sizeof(int));
    memcpy(this->nzval, nzval, n * M * sizeof(double));
}

SparseMatrix::SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzval, int m, int n)
{   
    assert(m > 0);
    assert(n > 0);

    this->nnz = nnz;    
    this->m = m;
    this->n = n;

    // construct matrix
    int i, j;
    double **mat = (double **) calloc(m, sizeof(double *));
    for (i = 0; i < m; i++)
        mat[i] = (double *) calloc(n, sizeof(double));

    for (i = 0; i < nnz; i++)
        mat[ridx[i]][cidx[i]] = nzval[i];
    
    // find M
    M = 0;
    for (j = 0; j < n; j++){
        for (i = 0; i < m; i++){
            M = (mat[i][j] != 0.0 && i > M) ? i : M;
        }
    }  

    // convert to ELL
    rowidx = new int[n * M];
    for (i = 0; i < n * M; i++)
        rowidx[i] = -1;
    this->nzval = new double[n * M];
    int nnzInCol;
    for (j = 0; j < n; j++){
        nnzInCol = 0;
        for (i = 0; i < m; i++){
            if (mat[i][j] != 0.0){
                rowidx[j * M + nnzInCol] = i;
                this->nzval[j * M + nnzInCol] = mat[i][j];
                nnzInCol++;
            }
        }
    }

    // free matrix
    for (int i = 0; i < m; i++)
        free(mat[i]);
    free(mat);
}

SparseMatrix::~SparseMatrix()
{
    delete[] rowidx;
    delete[] nzval;
}

int SparseMatrix::GetSize(int which) const
{
    assert(which == 1 || which == 2);
    if (which == 1)
        return m;
    return n;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix& otherMatrix)
{   
    if (n * M != otherMatrix.n * otherMatrix.M){
        delete[] rowidx;
        rowidx = new int[otherMatrix.n * otherMatrix.M];
        delete[] nzval;
        nzval = new double[otherMatrix.n * otherMatrix.M];
    }
    M = otherMatrix.M;
    m = otherMatrix.GetSize(1);
    n = otherMatrix.GetSize(2);
    nnz = otherMatrix.nnz;
    
    memcpy(rowidx, otherMatrix.rowidx, n * M * sizeof(int));
    memcpy(nzval, otherMatrix.nzval, n * M * sizeof(double));

    return *this;
}

SparseMatrix SparseMatrix::operator+() const
{
    SparseMatrix res(m, n, nnz, M, rowidx, nzval);
    return res;
}

SparseMatrix SparseMatrix::operator-() const
{
    SparseMatrix res(m, n, nnz, M, rowidx, nzval);
    for (int i = 0; i < n * M; i++)
        (res.nzval)[i] = (res.nzval)[i] * -1;
    return res;
}

SparseMatrix SparseMatrix::operator*(double a) const
{
    SparseMatrix res(m, n, nnz, M, rowidx, nzval);
    for (int i = 0; i < n * M; i++)
        (res.nzval)[i] = (res.nzval)[i] * a;
    return res;
}

Vector operator*(const SparseMatrix& m, const Vector& v)
{
    assert(m.n == v.GetSize());

    Vector res(m.m);
    int i, j, k;

    for (i = 0; i < m.n; i++){
        for (j = 0; j < m.M; j++){
            k = i * m.M + j;
            if (m.rowidx[k] == -1)
                break;
            res[m.rowidx[k]] += m.nzval[k] * v.Read(i);
        }
    }
    return res;
}
