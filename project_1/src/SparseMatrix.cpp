#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(int numRows, int numCols)
{
    m = numRows;
    n = numCols;
    M = 0;
    nnz = 0;
    rowidx = NULL;
    nzval = NULL;
}

SparseMatrix::SparseMatrix(int m, int n, int nnz, int M, int* rowidx, double* nzval)
{
    this->m = m;
    this->n = n;
    this->M = M;
    this->nnz = nnz;
    if (this->rowidx != NULL)
        delete[] this->rowidx;
    this->rowidx = new int[n * M];
    if (this->nzval != NULL)
        delete[] this->nzval;
    this->nzval = new double[n * M];
    for (int i = 0; i < n * M; i++){
        this->rowidx[i] = rowidx[i];
        this->nzval[i] = nzval[i];
    }
}

SparseMatrix::SparseMatrix(int nnz, int const *ridx, int const *cidx, double const *nzval, int m, int n)
{   
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
    if (rowidx != NULL)
        delete[] rowidx;
    if (nzval != NULL)
        delete[] nzval;
}

int SparseMatrix::GetSize(int which) const
{
    // assert(which == 1 || which == 2);
    if (which == 1)
        return m;
    return n;
}

SparseMatrix& SparseMatrix::operator=(const SparseMatrix& otherMatrix)
{   
    if (n * M != otherMatrix.n * otherMatrix.M){
        if (rowidx != NULL)
            delete[] rowidx;
        rowidx = new int[otherMatrix.n * otherMatrix.M];
        if (nzval != NULL)
            delete[] nzval;
        nzval = new double[otherMatrix.n * otherMatrix.M];
    }
    M = otherMatrix.M;
    m = otherMatrix.GetSize(1);
    n = otherMatrix.GetSize(2);
    nnz = otherMatrix.nnz;
    for (int i = 0; i < n * M; i++){
        rowidx[i] = otherMatrix.rowidx[i];
        nzval[i] = otherMatrix.nzval[i];
    }
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