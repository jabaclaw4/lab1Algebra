#include "../headers/matrix.h"
#include <cstdlib>
#include <iostream>

Matrix generate_random_matrix(int n, unsigned int seed) {
    srand(seed);
    Matrix A(n, Vector(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            //rand() даёт число от 0 до RAND_MAX делим на RAND_MAX получаем число от 0 до 1,
            // умножаем на 2 получаем от 0 до 2, вычитаем 1  получаем от -1 до 1 (раввномерное распределение от -1 до 1)
            A[i][j] = 2.0 * rand() / RAND_MAX - 1.0;
    return A;
}

Vector generate_random_vector(int n, unsigned int seed) {
    srand(seed);
    Vector b(n);
    for (int i = 0; i < n; i++)
        b[i] = 2.0 * rand() / RAND_MAX - 1.0;
    return b;
}

Matrix generate_hilbert(int n) {
    Matrix H(n, Vector(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            H[i][j] = 1.0 / (i + j + 1);
    return H;
}



void print_matrix(const Matrix& A) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            std::cout << A[i][j] << " ";
        std::cout << "\n";
    }
}

void print_vector(const Vector& v) {
    for (double x : v)
        std::cout << x << " ";
    std::cout << "\n";
}