#ifndef MATRIX_H
#define MATRIX_H

#include "Vector.hpp"
#include "Utilities.hpp"

typedef unsigned int uint;

template<typename T>
class Matrix{


private:
    Vector<T> matrix;
    uint rows;
    uint columns;
};



#endif