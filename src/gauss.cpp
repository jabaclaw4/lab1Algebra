#include "../headers/gauss.h"
#include <cmath>
#include <stdexcept>

//обратный ход для уже приведённой верхнетреугольной системы
static Vector back_substitution(const Matrix& A, const Vector& b) {
    int n = A.size();
    Vector x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = b[i];
        for (int j = i + 1; j < n; j++)
            x[i] -= A[i][j] * x[j];
        x[i] /= A[i][i];
    }
    return x;
}

//метод гаусса без выбора ведущего элемента
Vector gauss_no_pivot(Matrix A, Vector b) {
    int n = A.size();

    //прямой ход: приводим матрицу к верхнетреугольному виду
    for (int k = 0; k < n; k++) {
        if (fabs(A[k][k]) < 1e-15)
            throw std::runtime_error("zero pivot element");

        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++)
                A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    return back_substitution(A, b);
}

//метод гаусса с частичным выбором ведущего элемента по столбцу
Vector gauss_pivot(Matrix A, Vector b) {
    int n = A.size();

    for (int k = 0; k < n; k++) {
        //ищем строку с максимальным по модулю элементом в столбце k
        int max_row = k;
        for (int i = k + 1; i < n; i++)
            if (fabs(A[i][k]) > fabs(A[max_row][k]))
                max_row = i;

        //переставляем строки если нужно
        if (max_row != k) {
            std::swap(A[k], A[max_row]);
            std::swap(b[k], b[max_row]);
        }

        if (fabs(A[k][k]) < 1e-15)
            throw std::runtime_error("singular matrix");

        for (int i = k + 1; i < n; i++) {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++)
                A[i][j] -= factor * A[k][j];
            b[i] -= factor * b[k];
        }
    }

    return back_substitution(A, b);
}