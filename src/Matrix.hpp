#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.hpp"
#include "Utilities.hpp"

typedef unsigned int uint;

template<typename T>
class Matrix{
    Matrix(uint rows, uint cols, Vector<T> matrix)
        : rows(rows), cols(cols), matrix(matrix) 
    { }

    Matrix(uint size)
        : rows(size), cols(size) {
        Vector<T> identity(size*size, 0);
        for(int i=0; i<size; ++i)
            identity[i*(size+1)] = 1;
        matrix = identity;
    }

    inline uint rowNumber() const { return rows; }
    inline uint colNumber() const { return cols; }
    inline T get(uint row, uint column) const { return matrix[row*cols + column]; }
    inline T& get(uint row, uint column) { return matrix[row*cols + column]; }
    inline int index(uint row, uint column) const { return row*cols + column; }

private:
    Vector<T> matrix;
    uint rows;
    uint cols;
};



#endif