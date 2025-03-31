#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include "Matrix.hpp"
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

    id.print();
    
    Matrix<bool> Q(4);
    Matrix<bool> Qinv(4);

    Matrix<bool>::getRowEchelonForm(id, Q, Qinv);
    id.print();

    return 0;
}

