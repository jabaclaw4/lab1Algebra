#ifndef LU_H
#define LU_H

#include "matrix.h"

//структура для хранения результата LU-разложения (LU-разложение возвращает две матрицы L и U которые нужно сохранить отдельн чтобы
//потом использовать их много раз для разных правых частей)
//A = L * Uгде L нижняя треугольная с единичной диагональю, U  верхняя
struct LU {
    Matrix L;
    Matrix U;
};

//LU-разложение матрицы A без перестановок
LU lu_decompose(const Matrix& A);

//прямая подстановка: решает L*y = b
Vector forward_sub(const Matrix& L, const Vector& b);

//обратная подстановка: решает U*x = y
Vector back_sub(const Matrix& U, const Vector& y);

//решение системы через готовое LU-разложение
Vector lu_solve(const LU& lu, const Vector& b);

#endif