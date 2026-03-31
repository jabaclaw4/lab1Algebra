#ifndef UTILS_H
#define UTILS_H
#include "matrix.h"
#include "ctime"

//замер времени через clock()
struct Timer {
    clock_t start;
    void begin() {
        start = clock();
    }

    //возвращает время в миллисекундах
    double elapsed_ms() const {
        return 1000.0 * (clock() - start) / CLOCKS_PER_SEC;
    }
};

//евклидова норма вектора
double norm(const Vector& v);

//невязка: ||A*x - b||
double residual(const Matrix& A, const Vector& x, const Vector& b);

//относительная погрешность: ||x_approx - x_exact|| / ||x_exact||
double relative_error(const Vector& x_approx, const Vector& x_exact);

#endif