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

    inline uint rowNumber() const { return rows; }
    inline uint colNumber() const { return cols; }
    inline uint index(const uint row, const uint column) const { return row*cols + column; }

    inline const Vector<T>& data() const { return matrix; }
    inline T get(const uint row, uint const column) const { return matrix[row*cols + column]; }
    inline T& get(const uint row, const uint column) { return matrix[row*cols + column]; }
    Vector<T> getRow(const uint row) const;
    Vector<T> getColumn(const uint column) const;
    inline void nextRow(uint& index) const { index+=cols; }
    inline void nextColumn(uint& index) const { ++index; }
    void print(std::string endChar = "\n", std::ostream& outputStream = std::cout) const;
    Matrix slice(const uint rowBegin, const uint rowEnd, const uint columnBegin, const uint columnEnd) const;
    Matrix static transpose(const Matrix<T>& base);

    inline T operator()(const uint row, const uint column) const { return matrix[row*cols + column]; }
    inline T& operator()(const uint row, const uint column) { return matrix[row*cols + column]; }

    template<typename U>
    friend Matrix<U> operator+(const Matrix<U>& a, const Matrix<U>& b);

    void swapRows(uint rowA, uint rowB);
    void swapColumns(uint columnA, uint columnB);
    void changeRowSign(uint row);
    void changeColumnSign(uint column);
    void addRowMultiple(uint row, uint addedRow, T multiple);
    void addColumnMultiple(uint column, uint addedColumn, T multiple);

    void static swapRowsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint rowA, const uint rowB);
    void static swapColumnsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint columnA, const uint columnB);
    void static changeRowSignOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row);
    void static changeColumnSignOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint column);
    void static addRowMultipleOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row, const uint addedRow, const T multiple);
    void static addColumnMultipleOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint addedColumn, const uint column, const T multiple);

    void static reduceRowPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row, const uint column);
    void static reduceColumnPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Q_, const uint row, const uint column);

private:
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
void Matrix<T>::print(std::string endChar, std::ostream& out) const{
    for(uint i=0; i<rows; ++i){
        for(uint j=0; j<cols; ++j)
            out << (*this)(i, j) << ' ';
        out << '\n';
    }
    out << endChar;
}

template<typename T>
Matrix<T> Matrix<T>::slice(uint rowBegin, uint rowEnd, uint columnBegin, uint columnEnd) const{
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
    B.exchangeRows(rowA, rowB);
    Qinv.exchangeRows(rowA, rowB);
    Q.exchangeColumns(rowA, rowB);
}

template<typename T>
void Matrix<T>::swapColumnsOperation(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint columnA, const uint columnB){
    B.exchangeColumns(columnA, columnB);
    Q.exchangeColumns(columnA, columnB);    
    Qinv.exchangeRows(columnA, columnB);
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
void Matrix<T>::reduceRowPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column){
    T s = B.get(row, column);
    for(uint i=row+1; i<B.rowNumber(); ++i){
        T r = B.get(i, column);
        Matrix::addRowMultipleOperation(B, Q, Qinv, i, row, -floor(r, s));
    }
}

template<typename T>
void Matrix<T>::reduceColumnPartially(Matrix<T>& B, Matrix<T>& Q, Matrix<T>& Qinv, const uint row, const uint column){
    T s = B.get(row, column);
    for(uint i=column+1; i<B.colNumber(); ++i){
        T r = B.get(row, i);
        Matrix::addColumnMultipleOperation(B, Q, Qinv, column, i, -floor(r, s));
    }
}

#endif