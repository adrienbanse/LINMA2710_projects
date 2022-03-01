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

    if (rowidx == NULL) 
        this->rowidx = NULL;
    if (nzval == NULL)
        this->nzval = NULL;
    this->rowidx = new int[nnz];
    this->nzval = new double[nnz];
    for (int i = 0; i < nnz; i++){
        this->rowidx[i] = rowidx[i];
        this->nzval[i] = nzval[i];
    }
}

SparseVector::SparseVector() : AbstractVector(1) // TODO: check this
{
    this->nnz = 0;
    this->rowidx = NULL;
    this->nzval = NULL;
}

SparseVector::~SparseVector()
{   
    if (rowidx != NULL)
        delete[] rowidx;
    if (nzval != NULL)
        delete[] nzval;
}

double SparseVector::Read(int i) const
{    
    // assert(i > -1);
    // assert(i < GetSize());
    for (int j = 0; j < nnz; j++){
        if(rowidx[j] == i){
            return nzval[j];
        }
    }
    return 0.0;
}

SparseVector& SparseVector::operator=(const SparseVector& otherVector)
{   
    mSize = otherVector.GetSize();
    if (nnz != otherVector.nnz){
        if (rowidx != NULL)
            delete[] rowidx;
        if (nzval != NULL)
            delete[] nzval;
        nnz = otherVector.nnz;
        rowidx = new int[nnz];
        nzval = new double[nnz];
    }
    for(int i = 0; i < nnz; i++){
        rowidx[i] = otherVector.rowidx[i];
        nzval[i] = otherVector.nzval[i];
    }
    return *this;
}

SparseVector SparseVector::operator+(const SparseVector& v1) const
{
    // assert(GetSize() == v1.GetSize());
    int* newRowIdx = new int[GetSize()];
    double* newNzVal = new double[GetSize()];
    int newNnz = 0;

    int i = 0, j = 0; // two pointers
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
    for (int i = 0; i < newNnz; i++){
        finalRowIdx[i] = newRowIdx[i];
        finalNzVal[i] = newNzVal[i];
    }
    if (newRowIdx != NULL)
        delete[] newRowIdx;
    if (newNzVal != NULL)
    delete[] newNzVal;

    SparseVector res(newNnz, finalRowIdx, finalNzVal, GetSize());

    if (finalRowIdx != NULL)
        delete[] finalRowIdx;
    if (finalNzVal != NULL)
        delete[] finalNzVal; 

    return res;
}