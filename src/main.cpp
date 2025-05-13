#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
#include "CubicalComplex.hpp"
#include <vector>

using namespace std;

int main(){
    srand(time(NULL));

    const uint N = 5; 
    const uint M = 6;
    Matrix<int> id(N, M, Vector<int>(N*M, N*M, 0));

    for(int i=0; i<N; ++i){
        for(int j=0; j<M; ++j)
            id(i,j) = rand() % 10 -4;
    }

    Matrix<int> r(N);
    Matrix<int> rinv(N);
    Matrix<int> c(M);
    Matrix<int> cinv(M);

    id.print();
    // r.print();
    // cinv.print();

    SmithForm<int> smf = id.getSmithForm();
    smf.matrix.print();

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

