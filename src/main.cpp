#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
#include "CubicalComplex.hpp"
#include <vector>

using namespace std;

int main(){

    Cube skeletonCube("[0,1]x[0,1]");
    CubicalSet skeleton(Vector<Cube>{skeletonCube});
    CubicalComplexZ2 complex(skeleton);

    Vector<Matrix<bool>> bd = complex.getMatrixOfBoundaryOperator();
    for(const auto& bd_k : bd)
        bd_k.print();

    return 0;
}

