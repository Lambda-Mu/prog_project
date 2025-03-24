#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
#include <vector>

using namespace std;

int main(){

    Matrix<bool> id(7);

    id.print();

    id.addRowMultiple(3, 2, 5);

    Matrix<bool> id2(1);

    id2 = Matrix<bool>::transpose(id);

    Matrix<bool>::transpose(id).print();

    (id + id2).print();

    return 0;
}

