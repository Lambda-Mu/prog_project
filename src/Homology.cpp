#include "Homology.hpp"

using namespace std;

template<>
void HomologyGroup<int>::printFundamentalDecomposition() const{
    if(bettiNumber)
        cout << "Z^" << bettiNumber;
    for(int i=0; i<torsionCoefficients.size(); ++i)
        cout << " + Z_" << torsionCoefficients[i];
    if(bettiNumber == 0 and torsionCoefficients.empty())
        cout << "0";
    cout << endl;
}

template<>
void HomologyGroup<bool>::printFundamentalDecomposition() const{
    cout << "betti number: " << bettiNumber << endl;
    cout << "torsion: " << torsionCoefficients << endl;
}