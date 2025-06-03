#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
#include "CubicalComplex.hpp"
#include "Homology.hpp"
#include <vector>

using namespace std;

int main(){

    Cube skeletonCube("[0,1]x[0,1]");
    CubicalSet skeleton(Vector<Cube>{skeletonCube});
    CubicalComplexZ2 complex(skeleton);

    Vector<Matrix<bool>> bd = complex.getMatrixOfBoundaryOperator();

    for(uint i=0; i<=2; ++i){
        cout << complex.base(i) << endl;
    }
    cout << endl;

    for(const auto& bd_k : bd){
        cout << "Matrix:\n"; 
        bd_k.print();
    }

    auto H = getHomologyGroupOfChainComplex(bd);

    return 0;
}

