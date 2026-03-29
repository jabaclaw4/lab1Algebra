#include "../headers/lu.h"
#include <cmath>
#include <stdexcept>

//LU-разложение по алгоритму дулиттла: диагональ L состоит из единиц
LU lu_decompose(const Matrix& A) {
    int n = A.size();
    Matrix L(n, Vector(n, 0.0));
    Matrix U(n, Vector(n, 0.0));

    for (int i = 0; i < n; i++) {
        //заполняем строку U
        for (int j = i; j < n; j++) {
            U[i][j] = A[i][j];
            for (int k = 0; k < i; k++)
                U[i][j] -= L[i][k] * U[k][j];
        }

        if (fabs(U[i][i]) < 1e-15)
            throw std::runtime_error("zero diagonal element, LU without pivoting is not possible");

        //заполняем столбец L
        L[i][i] = 1.0;
        for (int j = i + 1; j < n; j++) {
            L[j][i] = A[j][i];
            for (int k = 0; k < i; k++)
                L[j][i] -= L[j][k] * U[k][i];
            L[j][i] /= U[i][i];
        }
    }

    return {L, U};
}

//прямая подстановка: L*y = b, диагональ L единичная
Vector forward_sub(const Matrix& L, const Vector& b) {
    int n = L.size();
    Vector y(n, 0.0);
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++)
            y[i] -= L[i][j] * y[j];
    }
    return y;
}

//обратная подстановка: U*x = y
Vector back_sub(const Matrix& U, const Vector& y) {
    int n = U.size();
    Vector x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= U[i][j] * x[j];
        x[i] /= U[i][i];
    }
    return x;
}

//полное решение: сначала L*y = b, потом U*x = y
Vector lu_solve(const LU& lu, const Vector& b) {
    Vector y = forward_sub(lu.L, b);
    return back_sub(lu.U, y);
}