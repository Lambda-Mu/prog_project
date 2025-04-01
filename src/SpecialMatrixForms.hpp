#ifndef SPECIAL_MATRIX_FORMS
#define SPECIAL_MATRIX_FORMS

#include "Matrix.hpp"

typedef unsigned int uint;

template<typename T>
class RowEchelon{
    RowEchelon(Matrix<T>& B)
        : B(B) {
            Q = Matrix<T>(B.rows);
            Qinv = Matrix<T>(B.rows);
            getRowEchelonForm(B, Q, Qinv, lastNonzeroRow);
        }

    const Matrix<T> B;
    const Matrix<T> Q;
    const Matrix<T> Qinv;
    const uint lastNonzeroRow;
};

template<typename T>
class ColumnEchelon{
    ColumnEchelon(Matrix<T>& B)
        : B(B) {
            R = Matrix<T>(B.rows);
            Rinv = Matrix<T>(B.rows);
            getColumnEchelonForm(B, R, Rinv, lastNonzeroColumn);
        }

private:
    const Matrix<T> B;
    const Matrix<T> R;
    const Matrix<T> Rinv;
    const uint lastNonzeroColumn;
};

#endif