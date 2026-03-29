#include "../headers/utils.h"
#include <cmath>

double norm(const Vector& v) {
    double s = 0.0;
    for (double x : v)
        s += x * x;
    return sqrt(s);
}

double residual(const Matrix& A, const Vector& x, const Vector& b) {
    int n = A.size();
    Vector r(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            r[i] += A[i][j] * x[j];
        r[i] -= b[i];
    }
    return norm(r);
}

double relative_error(const Vector& x_approx, const Vector& x_exact) {
    int n = x_approx.size();
    Vector diff(n);
    for (int i = 0; i < n; i++)
        diff[i] = x_approx[i] - x_exact[i];
    return norm(diff) / norm(x_exact);
}