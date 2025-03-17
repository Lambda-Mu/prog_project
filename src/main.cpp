#include <iostream>

#include "Vector.hpp"
#include "Utilities.hpp"
#include <vector>

int main(){
    Vector<int> y;

    cout << y[0] << endl << endl;
    
    for(int i=0; i<20; ++i){
        y.pushBack(i);
        cout << y.size() << ", " << y.capacity() << endl;
        for(int j=0; j<=i; ++j)
            cout << y[j] << ' ';
        cout << endl << endl;
    }

    for(int i: y)
        cout << y[i] << ' ';
    cout << endl;

    std::string k = "BLABLABLA";
    Vector<std::string> x(10, k);
    cout << &k << endl;
    cout << &x[1] << endl;
    cout << k << endl;
    cout << x[1] << endl;

    Vector<int> z{1,3};

    cout << endl;
    cout << z;
    cout << endl;

    for(int i=0; i<z.size(); ++i)
        cout << z[i] << " ";
    cout << endl;

    std::vector<int> p;
    p.resize(100);

    cout << p[50] << endl;

    return 0;
}

