#ifndef LEFSCHETZ_COMPLEX_H
#define LEFSCHETZ_COMPLEX_H

#include "Vector.hpp"
#include "Matrix.hpp"

template<class Cell, typename M>
class LefschetzComplex{
public:
    // virtual M operator()(const Cell& a, const Cell& b) const = 0;
    // virtual Vector<Cell> getBase(int dimension) const = 0;
    // virtual Vector<M> getCoordinates(const Vector<Cell>& cells) const = 0;
    virtual Vector<Matrix<M>> getMatrixOfBoundaryOperator() const = 0;
    // virtual void reduce(const Cell& cell, const Cell& face) = 0;
};



#endif
