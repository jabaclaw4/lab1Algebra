#ifndef MATRIX_H
#define MATRIX_H
#include <vector>

//тип для матрицы: вектор строк
using Matrix = std::vector<std::vector<double>>;

//тип для вектора-столбца
using Vector = std::vector<double>;

//генерация случайной матрицы n×n, элементы из [-1, 1]
Matrix generate_random_matrix(int n, unsigned int seed = 42);

//генерация случайного вектора длины n
Vector generate_random_vector(int n, unsigned int seed = 137);

//построение матрицы гильберта размера n×n: H[i][j] = 1/(i+j-1)
Matrix generate_hilbert(int n);

void print_matrix(const Matrix& A);
void print_vector(const Vector& v);

#endif