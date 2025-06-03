#ifndef HOMOLOGY_HPP
#define HOMOLOGY_HPP

#include <set>
#include "Matrix.hpp"
#include "LefschetzComplex.hpp"
#include "Vector.hpp"

template<typename T>
class HomologyGroup{
public:
    HomologyGroup(uint bettiNumber, Vector<T> torsionCoefficients) :
        bettiNumber(bettiNumber), torsionCoefficients(torsionCoefficients) {}
    const uint bettiNumber;
    const Vector<uint> torsionCoefficients;

    void printFundamentalDecomposition() const;
    bool isTrivial() const { return !bettiNumber and torsionCoefficients.empty(); }
};

template<typename T>
Vector<QuotientGroup<T>> getHomologyGroupOfChainComplex(const Vector<Matrix<T>>& boundaryOperator){
    Vector<KernelImage<T>> kerIm;
    kerIm.reserve(boundaryOperator.size());
    for(uint i=0; i<boundaryOperator.size(); ++i){
        kerIm.push_back(boundaryOperator[i].getKernelAndImage());
    }
    
    Vector<QuotientGroup<T>> homology;
    homology.reserve(boundaryOperator.size());
    
    uint limit = boundaryOperator.size() - 1;
    for(uint i=0; i<limit; ++i){    
        homology.push_back(getQuotientGroup(kerIm[i].kernel, kerIm[i+1].image));
    }
    uint rows = kerIm[limit].kernel.rowNumber();
    Matrix<T> zeroMatrix(rows, 1, Vector<T>(rows, rows, 0));
    homology.push_back(getQuotientGroup(kerIm[limit].kernel, zeroMatrix));
    
    return homology;
}

template<typename T, typename Cell>
Vector<HomologyGroup<T>> getGeneratorsOfHomology(
    const Vector<QuotientGroup<T>>& H)
{
    Vector<HomologyGroup<T>> homology;
    for(uint k=0; k<H.size(); ++k){
        if(H[k].quotientGenerator.isZero()){
            homology.pushBack(HomologyGroup(0, Vector<T>{}));
            continue;
        }
        uint s = H[k].numberOfZeroClasses;
        uint t = H[k].beginFiniteGenerators;
        uint m = H[k].inclusionMap.rowNumber();
        uint bettiNumber = 0;
        Vector<T> torsionCoeffs;
        for(uint j=s; j<m; ++j){
            if(j >= t)
                ++bettiNumber;
            else
                torsionCoeffs.pushBack(H[k].inclusionMap(j,j));   
        }
        homology.pushBack(HomologyGroup(bettiNumber, torsionCoeffs));
    }
}

template<typename T, typename Cell>
Vector<HomologyGroup<T>> getHomology(const LefschetzComplex<Cell, T> K){
    Vector<Matrix<T>> D = K.getMatrixOfBoundaryOperator();
    Vector<QuotientGroup<T>> H = getHomologyGroupOfChainComplex(D);
    Vector<HomologyGroup<T>> G = getGeneratorsOfHomology(H);
    return G;
}

#endif