#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
#include "CubicalComplex.hpp"
#include <vector>

using namespace std;

int main(){
    srand(time(NULL));

    const unsigned int N = 4; 
    Matrix<int> id(4);

    for(int i=0; i<4; ++i){
        for(int j=0; j<4; ++j)
            id(i,j) = rand() % 20 -4;
    }

    // Matrix<int> copy(id);
    // Matrix<int> rowBase(4);
    // Matrix<int> rowBaseInv(4);


    // cout << "matrix:\n";
    // id.print();
    // Matrix<int>::reduceRow(id, rowBase, rowBaseInv, 0, 0);
    // cout << "after transformation:\n";
    // id.print();
    // cout << "rowBase:\n";
    // rowBase.print();
    // cout << "rowBaseInv:\n";
    // rowBaseInv.print();
    // cout << "identity:" << (rowBase * rowBaseInv == Matrix<int>(4)) << '\n';
    // cout << "original:" << (rowBaseInv * id == copy) << '\n';
    // cout << "transformed:" << (rowBase * copy == id) << '\n';

//    id.print();
//    ColumnEchelon<int> ce(id);
//    ce.print();
//    RowEchelon<int> re(id);
//    re.print();
//    (ce * re).print();


    return 0;
}

