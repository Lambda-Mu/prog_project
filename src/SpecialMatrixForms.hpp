#ifndef SPECIAL_MATRIX_FORMS
#define SPECIAL_MATRIX_FORMS

#include "Matrix.hpp"

typedef unsigned int uint;

template<typename T>
class RowEchelon : public Matrix<T>{
public:
    RowEchelon(const Matrix<T>& M)
        {
            *this = M;
            Q = Matrix<T>(M.rowNumber());
            Qinv = Matrix<T>(M.rowNumber());
            Matrix<T>::getRowEchelonForm(*this, Q, Qinv, lastNonzeroRow);
        }

    Matrix<T> getOriginal() const{
        return Q * (*this);
    }

    Matrix<T> Q;
    Matrix<T> Qinv;
    uint lastNonzeroRow;
};

template<typename T>
class ColumnEchelon{
    ColumnEchelon(const Matrix<T>& M)
        : B(M), R(Matrix<T>(M.colNumber())), Rinv(Matrix<T>(M.colNumber())) {
            Matrix<T>::getColumnEchelonForm(B, R, Rinv, lastNonzeroColumn);
        }

private:
    const Matrix<T> B;
    const Matrix<T> R;
    const Matrix<T> Rinv;
    const uint lastNonzeroColumn;
};

#endif
