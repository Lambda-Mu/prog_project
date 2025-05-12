#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.hpp"
#include "Utilities.hpp"

typedef unsigned int uint;

template<typename T>
class Matrix{
public:
    Matrix() : rows(0), cols(0), matrix(Vector<T>()) { }

    Matrix(uint rows, uint cols, const Vector<T>& matrix)
        : rows(rows), cols(cols), matrix(matrix) 
    { }

    Matrix(uint size)
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
    inline uint index(uint row, uint column) const { return row*cols + column; }

    inline const Vector<T>& data() const { return matrix; }
    Vector<T> getRow(uint row) const;
    Vector<T> getColumn(uint column) const;
    inline void moveToNextRow(uint& index) const { index+=cols; }
    inline void moveToNextColumn(uint& index) const { ++index; }
    void print(const std::string endChar = "\n", std::ostream& outputStream = std::cout) const;
    Matrix slice(uint rowBegin, uint rowEnd, uint columnBegin, uint columnEnd) const;
    Matrix static transpose(const Matrix<T>& base);

    inline T operator()(uint row, uint column) const { return matrix[row*cols + column]; }
    inline T& operator()(uint row, uint column) { return matrix[row*cols + column]; }
    template<typename U>
    friend Matrix<U> operator+(const Matrix<U>& a, const Matrix<U>& b);
    template<typename U>
    friend Matrix<U> operator*(const Matrix<U>& a, const Matrix<U>& b);
    template<typename U>
    friend bool operator==(const Matrix<U>& a, const Matrix<U>& b);

    void swapRows(uint rowA, uint rowB);
    void swapColumns(uint columnA, uint columnB);
    void changeRowSign(uint row);
    void changeColumnSign(uint column);
    void addRowMultiple(uint row, uint addedRow, const T multiple);
    void addColumnMultiple(uint column, uint addedColumn, const T multiple);

    bool isRowZero(uint startRow, uint column) const;
    bool isColumnZero(uint startRow, uint column) const;

    void static getRowEchelonForm(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, uint& lastNonzeroRow);
    void static getColumnEchelonForm(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint& lastNonzeroColumn);

    void getKernelImage(Matrix<T>& kernel, Matrix<T>& image) const;

private:
    Vector<T> matrix;
    uint rows;
    uint cols;

    void static swapRowsOperation(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
        uint rowA, uint rowB);
    void static swapColumnsOperation(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
        uint columnA, uint columnB);
    void static changeRowSignOperation(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, uint row);
    void static changeColumnSignOperation(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint column);
    void static addRowMultipleOperation(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
        uint row, uint addedRow, const T multiple);
    void static addColumnMultipleOperation(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
        uint addedColumn, uint column, const T multiple);

    void static reduceRowValuesPartially(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv,
        uint row, uint startColumn);
    void static reduceColumnValuesPartially(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv,
        uint startRow, uint column);

    OneIndexedValue<T> findSmallestNonzeroInRow(uint row, uint startColumn) const;
    OneIndexedValue<T> findSmallestNonzeroInColumn(uint startRow, uint column) const;
    DoubleIndexedValue<T> findSmallestNonzeroInSubmatrix(uint firstDiagonalEntryIndex) const;

    void static prepareRow(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
        uint row, uint column);
    void static prepareColumn(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
        uint row, uint column);
    
    void static reduceRow(
        Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
        uint row, uint column);
    void static reduceColumn(
        Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
        uint row, uint column);

    void moveMinimalNonzero(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv,
        Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint diagonalEntryIndex);

    struct DivisibilityCheck{
    public:
        DivisibilityCheck(bool divisible, uint rowIndex, uint columnIndex, T divisor)
            : divisible(divisible), rowIndex(rowIndex), columnIndex(columnIndex), divisor(divisor)
        { }

        const bool divisible;
        uint rowIndex;
        uint columnIndex;
        const T divisor;
    };

    DivisibilityCheck checkForDivisibility(uint diagonalEntryIndex) const;

    void getPartialSmithForm(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv,
        Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint diagonalEntryIndex);

};

template<typename T>
class RowEchelon{
public:
    RowEchelon(const Matrix<T>& matrix)
        : matrix(matrix), lastNonzeroRow(0) {
        rowBase = Matrix<T>(matrix.rowNumber());
        rowBaseInv = Matrix<T>(matrix.rowNumber());
        Matrix<T>::getRowEchelonForm(matrix, rowBase, rowBaseInv, lastNonzeroRow);
    }

    Matrix<T> matrix;
    Matrix<T> rowBase;
    Matrix<T> rowBaseInv;
    uint lastNonzeroRow;
};

template<typename T>
class ColumnEchelon{
public:
    ColumnEchelon(const Matrix<T>& matrix)
        : matrix(matrix), lastNonzeroColumn(0) {
        columnBase = Matrix<T>(matrix.rowNumber());
        columnBaseInv = Matrix<T>(matrix.rowNumber());
        Matrix<T>::getColumnEchelonForm(matrix, columnBase, columnBaseInv, lastNonzeroColumn);
    }

    Matrix<T> matrix;
    Matrix<T> columnBase;
    Matrix<T> columnBaseInv;
    uint lastNonzeroColumn;
};

template<typename T>
class KernelImage{
public:
    KernelImage(const Matrix<T>& matrix)
        : kernel(Matrix<T>()), image(Matrix<T>()) {
        getKernelImage(matrix, kernel, image);
    }

    Matrix<T> kernel;
    Matrix<T> image;
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
    Vector<T> productData(a.rows * b.cols, a.rows * b.cols, 0);
    Matrix<T> product(a.rows, b.cols, productData);
    for(uint i=0; i<a.rows; ++i){
        for(uint j=0; j<b.cols; ++j){
            for(uint k=0; k<a.cols; ++k)
                product(i,j) += a(i,k) * b(k,j);
        }
    }
    return product;
}

template<typename T>
bool operator==(const Matrix<T>& a, const Matrix<T>& b){
    return a.rows == b.rows && a.cols == b.cols && a.matrix == b.matrix;
}

template<typename T>
Vector<T> Matrix<T>::getRow(uint row) const{
    return matrix.subVector(index(row,0), index(row,row+cols));
}

template<typename T>
Vector<T> Matrix<T>::getColumn(uint column) const{
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
void Matrix<T>::swapRows(uint rowA, uint rowB){
    for(uint i=0; i<cols; ++i){
        swap((*this)(rowA, i), (*this)(rowB, i));
    }
    return;
}

template<typename T>
void Matrix<T>::swapColumns(uint columnA, uint columnB){
    for(uint i=0; i<rows; ++i){
        swap((*this)(i, columnA), (*this)(i, columnB));
    }
    return;
}

template<typename T>
void Matrix<T>::changeRowSign(uint row){
    for(uint i=0; i<cols; ++i)
        (*this)(row, i) *= -1;
    return;
}

template<typename T>
void Matrix<T>::changeColumnSign(uint column){
    for(uint i=0; i<rows; ++i)
        (*this)(i, column) *= -1;
    return;
}

template<typename T>
void Matrix<T>::addRowMultiple(uint row, uint addedRow, const T multiple){
    for(uint i=0; i<cols; ++i)
        (*this)(row, i) += multiple * (*this)(addedRow, i);
}

template<typename T>
void Matrix<T>::addColumnMultiple(uint column, uint addedColumn, const T multiple){
    for(uint i=0; i<rows; ++i)
        (*this)(i, column) += multiple * (*this)(i, addedColumn);
}

template<typename T>
void Matrix<T>::swapRowsOperation(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
    uint rowA, uint rowB)
{
    matrix.swapRows(rowA, rowB);
    rowBase.swapRows(rowA, rowB);
    rowBaseInv.swapColumns(rowA, rowB);
}

template<typename T>
void Matrix<T>::swapColumnsOperation(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
    uint columnA, uint columnB)
{
    matrix.swapColumns(columnA, columnB);
    columnBase.swapColumns(columnA, columnB);    
    columnBaseInv.swapRows(columnA, columnB);
}

template<typename T>
void Matrix<T>::changeRowSignOperation(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
    uint row)
{
    matrix.changeRowSign(row);
    rowBase.changeRowSign(row);
    rowBaseInv.changeColumnSign(row);
}

template<typename T>
void Matrix<T>::changeColumnSignOperation(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
    uint column)
{
    matrix.changeColumnSign(column);
    columnBase.changeColumnSign(column);
    columnBaseInv.changeRowSign(column);
}

template<typename T>
void Matrix<T>::addRowMultipleOperation(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
    uint row, uint addedRow, const T multiple)
{
    matrix.addRowMultiple(row, addedRow, multiple);
    rowBase.addRowMultiple(row, addedRow, multiple);
    rowBaseInv.addColumnMultiple(addedRow, row, -multiple);
}

template<typename T>
void Matrix<T>::addColumnMultipleOperation(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
    uint column, uint addedColumn, const T multiple)
{
    matrix.addColumnMultiple(column, addedColumn, multiple);
    columnBase.addColumnMultiple(column, addedColumn, multiple);
    columnBaseInv.addRowMultiple(addedColumn, column, -multiple);
}

template<typename T>
void Matrix<T>::reduceRowValuesPartially(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
    uint row, uint startColumn)
{
    T s = matrix(row, startColumn);
    for(uint i=row+1; i<matrix.rowNumber(); ++i){
        T r = matrix(i, startColumn);
        Matrix::addRowMultipleOperation(matrix, rowBase, rowBaseInv, i, row, -floor(r, s));
    }
}

template<typename T>
void Matrix<T>::reduceColumnValuesPartially(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
    uint startRow, uint column)
{
    T s = matrix(startRow, column);
    for(uint i=column+1; i<matrix.colNumber(); ++i){
        T r = matrix(startRow, i);
        Matrix::addColumnMultipleOperation(matrix, columnBase, columnBaseInv, i, column, -floor(r, s));
    }
}

template<typename T>
OneIndexedValue<T> Matrix<T>::findSmallestNonzeroInRow(uint row, uint startColumn) const{
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
OneIndexedValue<T> Matrix<T>::findSmallestNonzeroInColumn(uint startRow, uint column) const{
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
DoubleIndexedValue<T> Matrix<T>::findSmallestNonzeroInSubmatrix(uint firstDiagonalEntryIndex) const{
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
bool Matrix<T>::isRowZero(uint row, uint startColumn) const{
    for(uint i=startColumn; i<cols; ++i){
        if((*this)(row, i) != 0)
            return false;
    }
    return true;
}

template<typename T>
bool Matrix<T>::isColumnZero(uint startRow, uint column) const{
    for(uint i=startRow; i<rows; ++i){
        if((*this)(i, column) != 0)
            return false;
    }
    return true;
}

template<typename T>
void Matrix<T>::prepareRow(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, 
    uint row, uint column)
{
    Matrix::swapRowsOperation(matrix, rowBase, rowBaseInv, 
        row, matrix.findSmallestNonzeroInColumn(row, column).index); 
}

template<typename T>
void Matrix<T>::prepareColumn(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
    uint row, uint column)
{
    Matrix::swapColumnsOperation(matrix, columnBase, columnBaseInv, 
        column, matrix.findSmallestNonzeroInRow(row, column).index); 
}

template<typename T>
void Matrix<T>::reduceRow(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv,
    uint row, uint column)
{
    while(!matrix.isColumnZero(row+1, column)){
        prepareRow(matrix, rowBase, rowBaseInv, row, column);
        reduceRowValuesPartially(matrix, rowBase, rowBaseInv, row, column);
    }
}

template<typename T>
void Matrix<T>::reduceColumn(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, 
    uint row, uint column)
{
    while(!matrix.isRowZero(row, column+1)){
        prepareColumn(matrix, columnBase, columnBaseInv, row, column);
        reduceColumnValuesPartially(matrix, columnBase, columnBaseInv, row, column);
    }
}

template<typename T>
void Matrix<T>::getRowEchelonForm(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv, uint& k){
    uint row = 0;
    uint column = 0;
    do{
        while(column < matrix.cols and matrix.isColumnZero(row, column))
            ++column;
        if(column == matrix.cols) 
            break;
        reduceRow(matrix, rowBase, rowBaseInv, row, column);
        ++row;
    }while(row < matrix.rows);
    k = row;
    return;
}

template<typename T>
void Matrix<T>::getColumnEchelonForm(Matrix<T>& matrix, Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint& k){
    int row = 0;
    int column = 0;
    do{
        while(row < matrix.rows and matrix.isRowZero(row, column)) 
            ++row;
        if(row == matrix.rows) 
            break;
        reduceColumn(matrix, columnBase, columnBaseInv, row, column);
        ++column;
    }while(column < matrix.cols);
    k = column;
    return;
}

template<typename T>
void Matrix<T>::getKernelImage(Matrix<T>& kernel, Matrix<T>& image) const{
    ColumnEchelon M(*this);
    uint k = M.lastNonzeroColumn;
    return KernelImage(M.columnBase.slice(0, cols, k+1, cols), M.matrix.slice(0, rows, 0, k));
}

template<typename T>
void Matrix<T>::moveMinimalNonzero(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv,
    Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint diagonalEntryIndex)
{
    DoubleIndexedValue position = matrix.findSmallestNonzeroInSubmatrix(diagonalEntryIndex);
    swapRowsOperation(matrix, rowBase, rowBaseInv, diagonalEntryIndex, position.rowIndex);
    swapColumnsOperation(matrix, columnBase, columnBaseInv, diagonalEntryIndex, position.columnIndex);
}

template<typename T>
typename Matrix<T>::DivisibilityCheck Matrix<T>::checkForDivisibility(uint diagonalEntryIndex) const{
    T value = (*this)(diagonalEntryIndex, diagonalEntryIndex);
    for(uint i=diagonalEntryIndex+1; i<rows; ++i){
        for(uint j=diagonalEntryIndex+1; j<cols; ++j){
            if((*this)(i,j) % value != 0)
                return DivisibilityCheck(false, i, j, floor((*this)(i,j), value));
        }
    }
    return DivisibilityCheck(true, 0, 0, 0);
}

template<typename T>
void Matrix<T>::getPartialSmithForm(Matrix<T>& matrix, Matrix<T>& rowBase, Matrix<T>& rowBaseInv,
    Matrix<T>& columnBase, Matrix<T>& columnBaseInv, uint diagonalEntryIndex)
{
    while(1){
        uint k = diagonalEntryIndex;
        moveMinimalNonzero(matrix, rowBase, rowBaseInv, columnBase, columnBaseInv, k);
        reducecolumnmatrixaseowPartially(matrix, rowBase, rowBaseInv, k, k);
        if(matrix.findSmallestNonzeroInColumn(k, k+1).value != 0)
            continue;
        reduceColumnPartially(matrix, columnBase, columnBaseInv, k, k);
        if(matrix.findSmallestNonzeroInColumn(k, k+1).value != 0)
            continue;
        DivisibilityCheck div = matrix.checkForDivisibility(k);
        if(div.divisible == 0){
            addRowMultipleOperation(matrix, rowBase, rowBaseInv, div.rowIndex, k, 1);
            addColumnMultipleOperation(matrix, columnBase, columnBaseInv, k, div.columnIndex, -div.quotient);
        }
        else return;
    }
}


#endif
