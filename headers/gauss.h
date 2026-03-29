#ifndef GAUSS_H
#define GAUSS_H

#include "matrix.h"

//метод гаусса без выбора ведущего элемента
Vector gauss_no_pivot(Matrix A, Vector b);

//метод гаусса с частичным выбором ведущего элемента по столбцу
Vector gauss_pivot(Matrix A, Vector b);

#endif