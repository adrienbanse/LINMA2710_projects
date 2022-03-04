#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "SparseVector.hpp"

SparseVector::SparseVector(int nnz, int const *rowidx, double const *nzval, int size) : AbstractVector(size)
{   
    this->nnz = nnz;
    this->rowidx = new int[nnz];
    this->nzval = new double[nnz];
    memcpy(this->rowidx, rowidx, nnz * sizeof(int));
    memcpy(this->nzval, nzval, nnz * sizeof(double));
}

SparseVector::SparseVector() : AbstractVector(1)
{
    this->nnz = 0;
    this->rowidx = new int[0];
    this->nzval = new double[0];
}

SparseVector::~SparseVector()
{   
    delete[] rowidx;
    delete[] nzval;
}

double SparseVector::Read(int i) const
{    
    assert(i > -1);
    assert(i < GetSize());
    for (int j = 0; j < nnz; j++){
        if(rowidx[j] == i){
            return nzval[j];
        }
    }
    return 0.;
}

SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{   
    mSize = otherVector.GetSize();
    if (nnz != otherVector.nnz){
        delete[] rowidx;
        delete[] nzval;
        nnz = otherVector.nnz;
        rowidx = new int[nnz];
        nzval = new double[nnz];
    }
    memcpy(rowidx, otherVector.rowidx, nnz * sizeof(int));
    memcpy(nzval, otherVector.nzval, nnz * sizeof(double));
    return *this;
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    assert(GetSize() == v1.GetSize());
    int* newRowIdx = new int[nnz + v1.nnz];
    double* newNzVal = new double[nnz + v1.nnz];
    int newNnz = 0;

    int i = 0, j = 0;
    while (i < nnz && j < v1.nnz){
        if (rowidx[i] < v1.rowidx[j]){
            newRowIdx[newNnz] = rowidx[i];
            newNzVal[newNnz] = nzval[i];
            i++;
            newNnz++;
        }
        else if (v1.rowidx[j] < rowidx[i]){
            newRowIdx[newNnz] = v1.rowidx[j];
            newNzVal[newNnz] = v1.nzval[j];
            j++;
            newNnz++;
        }
        else{
            newRowIdx[newNnz] = rowidx[i];
            newNzVal[newNnz] = nzval[i] + v1.nzval[j];
            i++;
            j++;
            newNnz++;
        }
    }
    while (i < nnz){
        newRowIdx[newNnz] = rowidx[i];
        newNzVal[newNnz] = nzval[i];
        i++;
        newNnz++;
    }
    while (j < v1.nnz){
        newRowIdx[newNnz] = v1.rowidx[j];
        newNzVal[newNnz] = v1.nzval[j];
        j++;
        newNnz++;
    }

    int* finalRowIdx = new int[newNnz];
    double* finalNzVal = new double[newNnz];
    memcpy(finalRowIdx, newRowIdx, newNnz * sizeof(int));
    memcpy(finalNzVal, newNzVal, newNnz * sizeof(double));

    delete[] newRowIdx;
    delete[] newNzVal; 

    SparseVector res(newNnz, finalRowIdx, finalNzVal, GetSize());

    delete[] finalRowIdx;
    delete[] finalNzVal; 

    return res;
}