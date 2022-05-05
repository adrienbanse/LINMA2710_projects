#include "SubMatrix.hpp"
#include <cassert>

double* SubMatrix::ptr(int i, int j){
    assert(i >= 1 && i <= numRows);
    assert(j >= 1 && j <= numCols);
    return &(mData[(i-1) * numCols + (j-1)]);
}

SubMatrix::SubMatrix(int numRows, int numCols){
    assert(numRows > 0);
    assert(numCols > 0);
    this->numRows = numRows;
    this->numCols = numCols;
    mData = new double[numRows * numCols];
    for (int i = 0; i < numRows * numCols; i ++)
        mData[i] = 0.;
    MPI_Type_vector(numRows, 1, numCols, MPI_DOUBLE, &stride);
    MPI_Type_commit(&stride);
}

SubMatrix::~SubMatrix() {
    delete[] mData;
    MPI_Type_free(&stride);
}

int SubMatrix::Size(int row) const{
    assert(row == 1 || row == 2);
    if (row == 1)
        return numRows;
    return numCols;
}

double SubMatrix::Read(int i, int j) const{
    assert(i >= 1 && i <= numRows);
    assert(j >= 1 && j <= numCols);
    return mData[(i-1) * numCols + (j-1)];
}

double& SubMatrix::operator()(int i, int j){
    assert(i >= 1 && i <= numRows);
    assert(j >= 1 && j <= numCols);
    return mData[(i-1) * numCols + (j-1)];
}

void SubMatrix::SendReceiveRows(int rowsend, int ranksend, int rowrecv, int rankrecv, int tag, MPI_Comm comm){
    MPI_Sendrecv(
        ptr(rowsend, 1), numCols, MPI_DOUBLE, ranksend, tag, 
        ptr(rowrecv, 1), numCols, MPI_DOUBLE, rankrecv, tag,
        comm, MPI_STATUS_IGNORE
    );
}

void SubMatrix::SendReceiveColumns(int colsend, int ranksend, int colrecv, int rankrecv, int tag, MPI_Comm comm){
    MPI_Sendrecv(
        ptr(1, colsend), 1, stride, ranksend, tag, 
        ptr(1, colrecv), 1, stride, rankrecv, tag,
        comm, MPI_STATUS_IGNORE
    );
}
