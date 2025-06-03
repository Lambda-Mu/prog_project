#include "Matrix.hpp"

template<>
Matrix<bool> operator+(const Matrix<bool>& a, const Matrix<bool>& b){
    uint size = a.rows * a.cols;
    Vector<bool> sum(size, size);
    for(uint i=0; i<size; ++i)
        sum[i] = a.matrix[i] xor b.matrix[i];
    return Matrix(a.rows, a.cols, sum);
}

template<>
Matrix<bool> operator*(const Matrix<bool>& a, const Matrix<bool>& b){
    Vector<bool> product(a.rows * b.cols, a.rows * b.cols, 0);
    for(uint i=0; i<a.rows; ++i){
        for(uint j=0; j<b.cols; ++j){
            for(uint k=0; k<a.cols; ++k)
                product[i * b.cols + j] = product[i * b.cols + j] xor a(i,k) xor b(k, j);
        }
    }
    return Matrix(a.rows, b.cols, product);
}

template<>
void Matrix<bool>::addRowMultiple(uint row, uint addedRow, const bool multiple){
    for(uint i=0; i<cols; ++i)
        (*this)(row, i) = (*this)(row, i) xor (*this)(addedRow, i);
}

template<>
void Matrix<bool>::addColumnMultiple(uint column, uint addedColumn, const bool multiple){
    for(uint i=0; i<rows; ++i)
        (*this)(i, column) = (*this)(i, column) xor (*this)(i, addedColumn);
}

template<>
void Matrix<bool>::reduceRowValuesPartially(Matrix<bool>& B, Matrix<bool>& Q, Matrix<bool>& Qinv, uint row, uint startColumn){
    bool s = B(row, startColumn);
    for(uint i=row+1; i<B.rowNumber(); ++i){
        if(B(i, startColumn))
            addRowMultipleOperation(B, Q, Qinv, i, row, 1);
    }
}

template<>
void Matrix<bool>::reduceColumnValuesPartially(Matrix<bool>& B, Matrix<bool>& Q, Matrix<bool>& Qinv, uint startRow, uint column){
    bool s = B(startRow, column);
    for(uint i=column+1; i<B.colNumber(); ++i){
        if(B(startRow, i))
            addColumnMultipleOperation(B, Q, Qinv, column, i, 1);
    }
}

template<>
OneIndexedValue<bool> Matrix<bool>::findSmallestNonzeroInRow(uint row, uint startColumn) const{
    for(uint k = startColumn; k<cols; ++k){
        if((*this)(row,k))
            return OneIndexedValue<bool>(k, 1);
    }
    return OneIndexedValue<bool>(startColumn, 0);
}

template<>
OneIndexedValue<bool> Matrix<bool>::findSmallestNonzeroInColumn(uint startRow, uint column) const{
    for(uint k = startRow; k<rows; ++k){
        if((*this)(k,column))
            return OneIndexedValue<bool>(k, 1);
    }
    return OneIndexedValue<bool>(startRow, 0);
}

template<>
DoubleIndexedValue<bool> Matrix<bool>::findSmallestNonzeroInSubmatrix(uint firstDiagonalEntryIndex) const{
    int value = (*this)(firstDiagonalEntryIndex, firstDiagonalEntryIndex);
    int rowIndex = firstDiagonalEntryIndex;
    int columnIndex = firstDiagonalEntryIndex;
    int dataIndex{};
    for(int i=firstDiagonalEntryIndex; i<rows; ++i){
        dataIndex = index(i, firstDiagonalEntryIndex);
        for(int j=firstDiagonalEntryIndex; j<cols; ++j){
            if( matrix[dataIndex] )
                return DoubleIndexedValue<bool>(1, i, j);
            ++dataIndex;
        }
    }
    return DoubleIndexedValue<bool>(value, rowIndex, columnIndex);
}

template<>
void Matrix<bool>::getPartialSmithForm(Matrix<bool>& matrix, Matrix<bool>& rowBase, Matrix<bool>& rowBaseInv,
    Matrix<bool>& columnBase, Matrix<bool>& columnBaseInv, uint diagonalEntryIndex)
{
    while(1){
        uint k = diagonalEntryIndex;
        moveMinimalNonzero(matrix, rowBase, rowBaseInv, columnBase, columnBaseInv, k);
        reduceRowValuesPartially(matrix, rowBase, rowBaseInv, k, k);
        if(matrix.findSmallestNonzeroInColumn(k, k+1).value != 0)
            continue;
        reduceColumnValuesPartially(matrix, columnBase, columnBaseInv, k, k);
        if(matrix.findSmallestNonzeroInRow(k, k+1).value != 0)
            continue;
        return;
    }
}

template<>
SmithForm<bool> Matrix<bool>::getSmithForm() const{
    Matrix<bool> B(*this);
    Matrix<bool> Q(rows);
    Matrix<bool> Qinv(rows);
    Matrix<bool> R(cols);
    Matrix<bool> Rinv(cols);
    if(cols == 1 && rows == 1){
        if((*this)(0,0) == 0)
            return SmithForm<bool>(B, Q, Qinv, R, Rinv, 0, 0);
        if((*this)(0,0) == 1)
            return SmithForm<bool>(B, Q, Qinv, R, Rinv, 1, 1);
        return SmithForm<bool>(B, Q, Qinv, R, Rinv, 0, 1);
    }
    uint numberOfOnes = 0;
    uint m = std::min(cols, rows);
    while(numberOfOnes < m and B.findSmallestNonzeroInSubmatrix(numberOfOnes).value != 0){
        getPartialSmithForm(B, Q, Qinv, R, Rinv, numberOfOnes);
        ++numberOfOnes;
    }
    return SmithForm(B, Q, Qinv, R, Rinv, numberOfOnes, numberOfOnes);
}