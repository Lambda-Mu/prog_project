#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.hpp"
#include "Utilities.hpp"

typedef unsigned int uint;

template<typename T>
class Matrix{
public:
    Matrix(const uint rows, const uint cols, const Vector<T>& matrix)
        : rows(rows), cols(cols), matrix(matrix) 
    { }

    Matrix(const uint size)
        : rows(size), cols(size) {
        Vector<T> identity(size*size, size*size, 0);
        for(uint i=0; i<size; ++i)
            identity[i*(size+1)] = 1;
        matrix = identity;
    }

    Matrix(const Matrix<T>& M)
        : rows(M.rows), cols(M.cols) {
            matrix = M.matrix;
        }

    inline uint rowNumber() const { return rows; }
    inline uint colNumber() const { return cols; }
    inline uint index(const uint row, const uint column) const { return row*cols + column; }

    inline const Vector<T>& data() const { return matrix; }
    Vector<T> getRow(const uint row) const;
    Vector<T> getColumn(const uint column) const;
    inline void moveToNextRow(uint& index) const { index+=cols; }
    inline void moveToNextColumn(uint& index) const { ++index; }
    void print(const std::string endChar = "\n", std::ostream& outputStream = std::cout) const;
    Matrix slice(const uint rowBegin, const uint rowEnd, const uint columnBegin, const uint columnEnd) const;
    Matrix static transpose(const Matrix<T>& base);

    inline T operator()(const uint row, const uint column) const { return matrix[row*cols + column]; }
    inline T& operator()(const uint row, const uint column) { return matrix[row*cols + column]; }
    template<typename U>
    friend Matrix<U> operator+(const Matrix<U>& a, const Matrix<U>& b);
    template<typename U>
    friend Matrix<U> operator*(const Matrix<U>& a, const Matrix<U>& b);

    void swapRows(const uint rowA, const uint rowB);
    void swapColumns(const uint columnA, const uint columnB);
    void changeRowSign(const uint row);
    void changeColumnSign(const uint column);
    void addRowMultiple(const uint row, const uint addedRow, const T multiple);
    void addColumnMultiple(const uint column, const uint addedColumn, const T multiple);

    void static swapRowsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint rowA, const uint rowB);
    void static swapColumnsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint columnA, const uint columnB);
    void static changeRowSignOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row);
    void static changeColumnSignOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint column);
    void static addRowMultipleOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row, const uint addedRow, const T multiple);
    void static addColumnMultipleOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint addedColumn, const uint column, const T multiple);

    void static reduceColumnValuesPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint startRow, const uint column);
    void static reduceRowValuesPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row, const uint startColumn);

    OneIndexedValue<T> findSmallestNonzeroInRow(const uint row, const uint startColumn) const;
    OneIndexedValue<T> findSmallestNonzeroInColumn(const uint startRow, const uint column) const;
    DoubleIndexedValue<T> findSmallestNonzeroInSubmatrix(const uint firstDiagonalEntryIndex) const;

    bool isRowZero(const uint startRow, const uint column) const;
    bool isColumnZero(const uint startRow, const uint column) const;

    void static prepareRow(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column);
    void static prepareColumn(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column);
    
    void static reduceRow(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column);
    void static reduceColumn(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column);

    void static getRowEchelonForm(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, uint& k);
    void static getColumnEchelonForm(Matrix<T>& B, Matrix<T>& R, Matrix<T>& Rinv, uint& k);

protected:
    Vector<T> matrix;
    uint rows;
    uint cols;
};

template<typename T>
Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b){
    uint size = a.rows * a.cols;
    Vector<T> sum(size, size);
    for(uint i=0; i<size; ++i)
        sum[i] = a.matrix[i] + b.matrix[i];
    return Matrix(a.rows, a.cols, sum);
}

template<typename T>
Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b){
    Vector<T> product(a.rows * b.cols, a.rows * b.cols, 0);
    for(uint i=0; i<a.rows; ++i){
        for(uint j=0; j<b.cols; ++j){
            for(uint k=0; k<a.cols; ++k)
                product[i * b.cols + j] += a(i,k) + b(k, j);
        }
    }
    return Matrix(a.rows, b.cols, product);
}

template<typename T>
Vector<T> Matrix<T>::getRow(const uint row) const{
    return matrix.subVector(index(row,0), index(row,row+cols));
}

template<typename T>
Vector<T> Matrix<T>::getColumn(const uint column) const{
    Vector<T> columnVec(rows);
    for(uint i=0; i<rows; ++i)
        columnVec.pushBack((*this)(i, column));
    return columnVec;
}

template<typename T>
void Matrix<T>::print(const std::string endChar, std::ostream& out) const{
    for(uint i=0; i<rows; ++i){
        for(uint j=0; j<cols; ++j)
            out << (*this)(i, j) << ' ';
        out << '\n';
    }
    out << endChar;
}

template<typename T>
Matrix<T> Matrix<T>::slice(const uint rowBegin, const uint rowEnd, const uint columnBegin, const uint columnEnd) const{
    Vector<T> submatrix;
    // if(columnEnd < columnBegin && rowEnd < rowBegin)
        // return Matrix(1, 1, Vector<T>{0});
    // if(columnEnd < columnBegin)
        // return Matrix(rowEnd-rowBegin+1, 1, Vector<T>(rowEnd-rowBegin+1, 0));
    // if(rowEnd < rowBegin)
        // return Matrix(1, columnEnd-columnBegin+1, Vector<T>(columnEnd-columnBegin+1, 0));

    submatrix.reserve((rowEnd-rowBegin)*(columnEnd-columnBegin));
    for(uint i=rowBegin; i<rowEnd; ++i)
        submatrix.copyBack(matrix, index(i, columnBegin), index(i, columnEnd));

    return Matrix(rowEnd-rowBegin, columnEnd-columnBegin, submatrix);
}

template<typename T>
Matrix<T> Matrix<T>::transpose(const Matrix<T>& base){
    Vector<T> transposed(base.cols*base.rows);
    for(uint i=0; i<base.cols; ++i){
        for(uint j=0; j<base.rows; ++j)
            transposed.pushBack(base(j,i));
    }

    return Matrix(base.cols, base.rows, transposed);
}

template<typename T>
void Matrix<T>::swapRows(const uint rowA, const uint rowB){
    for(uint i=0; i<cols; ++i){
        swap((*this)(rowA, i), (*this)(rowB, i));
    }
    return;
}

template<typename T>
void Matrix<T>::swapColumns(const uint columnA, const uint columnB){
    for(uint i=0; i<rows; ++i){
        swap((*this)(i, columnA), (*this)(i, columnB));
    }
    return;
}

template<typename T>
void Matrix<T>::changeRowSign(const uint row){
    for(uint i=0; i<cols; ++i)
        (*this)(row, i) *= -1;
    return;
}

template<typename T>
void Matrix<T>::changeColumnSign(const uint column){
    for(uint i=0; i<rows; ++i)
        (*this)(i, column) *= -1;
    return;
}

template<typename T>
void Matrix<T>::addRowMultiple(const uint row, const uint addedRow, const T multiple){
    for(uint i=0; i<cols; ++i)
        (*this)(row, i) += multiple * (*this)(addedRow, i);
}

template<typename T>
void Matrix<T>::addColumnMultiple(const uint column, const uint addedColumn, const T multiple){
    for(uint i=0; i<rows; ++i)
        (*this)(i, column) += multiple * (*this)(i, addedColumn);
}

template<typename T>
void Matrix<T>::swapRowsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint rowA, const uint rowB){
    B.swapRows(rowA, rowB);
    Qinv.swapRows(rowA, rowB);
    Q.swapColumns(rowA, rowB);
}

template<typename T>
void Matrix<T>::swapColumnsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint columnA, const uint columnB){
    B.swapColumns(columnA, columnB);
    Q.swapColumns(columnA, columnB);    
    Qinv.swapRows(columnA, columnB);
}

template<typename T>
void Matrix<T>::changeRowSignOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row){
    B.changeRowSign(row);
    Qinv.changeRowSign(row);
    Q.changeColumnSign(row);
}

template<typename T>
void Matrix<T>::changeColumnSignOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint column){
    B.changeColumnSign(column);
    Q.changeColumnSign(column);
    Qinv.changeRowSign(column);
}

template<typename T>
void Matrix<T>::addRowMultipleOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint addedRow, const T multiple){
    B.addRowMultiple(row, addedRow, multiple);
    Q.addColumnMultiple(row, addedRow, -multiple);
    Qinv.addRowMultiple(row, addedRow, multiple);
}

template<typename T>
void Matrix<T>::addColumnMultipleOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint column, const uint addedColumn, const T multiple){
    B.addColumnMultiple(column, addedColumn, multiple);
    Q.addColumnMultiple(column, addedColumn, multiple);
    Qinv.addRowMultiple(column, addedColumn, -multiple);
}

template<typename T>
void Matrix<T>::reduceRowValuesPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint startColumn){
    T s = B(row, startColumn);
    for(uint i=row+1; i<B.rowNumber(); ++i){
        T r = B(i, startColumn);
        Matrix::addRowMultipleOperation(B, Q, Qinv, i, row, -floor(r, s));
    }
}

template<typename T>
void Matrix<T>::reduceColumnValuesPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint startRow, const uint column){
    T s = B(startRow, column);
    for(uint i=column+1; i<B.colNumber(); ++i){
        T r = B(startRow, i);
        Matrix::addColumnMultipleOperation(B, Q, Qinv, column, i, -floor(r, s));
    }
}

template<typename T>
OneIndexedValue<T> Matrix<T>::findSmallestNonzeroInRow(const uint row, const uint startColumn) const{
    uint indexMin = startColumn;
    T value = absoluteValue((*this)(row, startColumn));
    for(uint i=startColumn+1; i<cols; ++i){
        if((absoluteValue((*this)(row, i)) < value and (*this)(row, i) != 0) or value == 0){
            value = absoluteValue((*this)(row, i));
            indexMin = i;
            if(value == 1)
                return OneIndexedValue<T>(indexMin, value);
        }
    }
    return OneIndexedValue<T>(indexMin, value);
}

template<typename T>
OneIndexedValue<T> Matrix<T>::findSmallestNonzeroInColumn(const uint startRow, const uint column) const{
    uint indexMin = startRow;
    T value = absoluteValue((*this)(startRow, column));
    for(uint i=startRow+1; i<rows; ++i){
        if((absoluteValue((*this)(i, column)) < value and (*this)(i, column) != 0) or value == 0){
            value = absoluteValue((*this)(i, column));
            indexMin = i;
        }
    }
    return OneIndexedValue<T>(indexMin, value);
}

template<typename T>
DoubleIndexedValue<T> Matrix<T>::findSmallestNonzeroInSubmatrix(const uint firstDiagonalEntryIndex) const{
    T value = absoluteValue((*this)(firstDiagonalEntryIndex, firstDiagonalEntryIndex));
    uint rowIndex = firstDiagonalEntryIndex;
    uint columnIndex = firstDiagonalEntryIndex;
    uint dataIndex{};
    for(int i=firstDiagonalEntryIndex; i<rows; ++i){
        dataIndex = index(i, firstDiagonalEntryIndex);
        for(int j=firstDiagonalEntryIndex; j<cols; ++j){
            if( (absoluteValue(matrix[dataIndex]) < value and matrix[dataIndex] != 0) or value == 0 ){
                value = matrix[dataIndex];
                rowIndex = i;
                columnIndex = j;
            }
            ++dataIndex;
        }
    }    
    return DoubleIndexedValue<T>(value, rowIndex, columnIndex);
}

template<typename T>
bool Matrix<T>::isRowZero(const uint row, const uint startColumn) const{
    for(uint i=startColumn; i<cols; ++i){
        if((*this)(row, i) != 0)
            return false;
    }
    return true;
}

template<typename T>
bool Matrix<T>::isColumnZero(const uint startRow, const uint column) const{
    for(uint i=startRow; i<rows; ++i){
        if((*this)(i, column) != 0)
            return false;
    }
    return true;
}

template<typename T>
void Matrix<T>::prepareRow(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column){
    Matrix::swapRowsOperation(B, Q, Qinv, row, B.findSmallestNonzeroInColumn(row, column).index); 
}

template<typename T>
void Matrix<T>::prepareColumn(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column){
    Matrix::swapColumnsOperation(B, Q, Qinv, column, B.findSmallestNonzeroInRow(row, column).index); 
}

template<typename T>
void Matrix<T>::reduceRow(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column){
    while(!B.isColumnZero(row+1, column)){
        prepareRow(B, Q, Qinv, row, column);
        reduceRowValuesPartially(B, Q, Qinv, row, column);
    }
}

template<typename T>
void Matrix<T>::reduceColumn(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column){
    while(!B.isRowZero(row, column+1)){
        prepareColumn(B, Q, Qinv, row, column);
        reduceColumnValuesPartially(B, Q, Qinv, row, column);
    }
}

template<typename T>
void Matrix<T>::getRowEchelonForm(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, uint& k){
    uint row = 0;
    uint column = 0;
    do{
        while(column < B.cols and B.isColumnZero(row, column))
            ++column;
        if(column == B.cols) 
            break;
        reduceRow(B, Q, Qinv, row, column);
        ++row;
    }while(row < B.rows);
    k = row;
    return;
}

template<typename T>
void Matrix<T>::getColumnEchelonForm(Matrix<T>& B, Matrix<T>& R, Matrix<T>& Rinv, uint& k){
    int row = 0;
    int column = 0;
    do{
        while(row < B.rows and B.isRowZero(row, column)) 
            ++row;
        if(row == B.rows) 
            break;
        reduceColumn(B, R, Rinv, row, column);
        ++column;
    }while(column < B.cols);
    k = column;
    return;
}

#endif