#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
#include "SpecialMatrixForms.hpp"
#include "CubicalComplex.hpp"
#include <vector>

using namespace std;

int main(){
    srand(time(NULL));

    const unsigned int N = 4; 
    Matrix<bool> id(4);

    for(int i=0; i<4; ++i){
        for(int j=0; j<4; ++j)
            id(i,j) = rand() %2;
    }

    // id.print();
    
    Matrix<bool> Q(4);
    Matrix<bool> Qinv(4);

    // Matrix<bool>::getRowEchelonForm(id, Q, Qinv);
    // id.print();

    // Matrix<int> k(N);
    // for(int i=0; i<N; ++i){
    //     for(int j=0; j<N; ++j)
    //         k(i,j) = ( rand() % 30 ) - 10;
    // }

    // k.print();

    // RowEchelon re(k);

    // re.print();

    // re.Q.print();
    // re.Qinv.print();

    // cout << re.lastNonzeroRow << endl;

    // re.getOriginal().print();

    // re.B.print();
    // cout << re.lastNonzeroRow << endl;

    return 0;
}

