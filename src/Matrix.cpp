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
void Matrix<bool>::addRowMultiple(const uint row, const uint addedRow, const bool multiple){
    for(uint i=0; i<cols; ++i)
        (*this)(row, i) = (*this)(row, i) xor (*this)(addedRow, i);
}

template<>
void Matrix<bool>::addColumnMultiple(const uint column, const uint addedColumn, const bool multiple){
    for(uint i=0; i<rows; ++i)
        (*this)(i, column) = (*this)(i, column) xor (*this)(i, addedColumn);
}

template<>
void Matrix<bool>::reduceRowValuesPartially(Matrix<bool>& B, Matrix<bool>& Q, Matrix<bool>& Qinv, const uint row, const uint startColumn){
    bool s = B(row, startColumn);
    for(uint i=row+1; i<B.rowNumber(); ++i){
        if(B(i, startColumn))
            addRowMultipleOperation(B, Q, Qinv, i, row, 1);
    }
}

template<>
void Matrix<bool>::reduceColumnValuesPartially(Matrix<bool>& B, Matrix<bool>& Q, Matrix<bool>& Qinv, const uint startRow, const uint column){
    bool s = B(startRow, column);
    for(uint i=column+1; i<B.colNumber(); ++i){
        if(B(startRow, i))
            addColumnMultipleOperation(B, Q, Qinv, column, i, 1);
    }
}

template<>
OneIndexedValue<bool> Matrix<bool>::findSmallestNonzeroInRow(const uint row, const uint startColumn) const{
    for(uint k = startColumn; k<cols; ++k){
        if((*this)(row,k))
            return OneIndexedValue<bool>(k, 1);
    }
    return OneIndexedValue<bool>(startColumn, 0);
}

template<>
OneIndexedValue<bool> Matrix<bool>::findSmallestNonzeroInColumn(const uint startRow, const uint column) const{
    for(uint k = startRow; k<rows; ++k){
        if((*this)(k,column))
            return OneIndexedValue<bool>(k, 1);
    }
    return OneIndexedValue<bool>(startRow, 0);
}

template<>
DoubleIndexedValue<bool> Matrix<bool>::findSmallestNonzeroInSubmatrix(const uint firstDiagonalEntryIndex) const{
    int value = (firstDiagonalEntryIndex, firstDiagonalEntryIndex);
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