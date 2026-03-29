#include <iostream>
#include <vector>
#include <ctime>
#include "headers/matrix.h"
#include "headers/gauss.h"
#include "headers/lu.h"
#include "headers/utils.h"

//эксперимент 1: сравнение времени для разных размеров матрицы
void experiment_1() {
    std::cout << "\n=== Experiment 1: solving time for different matrix sizes ===\n";
    std::cout << "n\tGauss (ms)\tGauss+pivot (ms)\tLU decomp (ms)\tLU solve (ms)\tLU total (ms)\n";
    std::cout << "-------------------------------------------------------------------------------------\n";

    std::vector<int> sizes = {100, 200, 500, 1000};

    for (int n : sizes) {
        Matrix A = generate_random_matrix(n, 40);
        Vector b = generate_random_vector(n, 41);

        Timer t;

        t.begin();
        gauss_no_pivot(A, b);
        double t_gauss = t.elapsed_ms();

        t.begin();
        gauss_pivot(A, b);
        double t_gauss_piv = t.elapsed_ms();

        double t_lu_decomp = -1, t_lu_solve = -1;
        try {
            t.begin();
            LU lu = lu_decompose(A);
            t_lu_decomp = t.elapsed_ms();

            t.begin();
            lu_solve(lu, b);
            t_lu_solve = t.elapsed_ms();
        } catch (...) {
            //матрица не допускает LU без перестановок
        }

        std::cout << n << "\t" << t_gauss << "\t\t" << t_gauss_piv << "\t\t\t";
        if (t_lu_decomp < 0)
            std::cout << "n/a\t\tn/a\t\tn/a\n";
        else
            std::cout << t_lu_decomp << "\t\t" << t_lu_solve << "\t\t" << t_lu_decomp + t_lu_solve << "\n";
    }
}

//эксперимент 2: несколько правых частей — преимущество LU
void experiment_2() {
    std::cout << "\n=== Experiment 2: multiple right-hand sides, n = 500 ===\n";
    std::cout << "k\tGauss+pivot (ms)\tLU (ms)\n";
    std::cout << "----------------------------------------\n";

    int n = 500;
    Matrix A = generate_random_matrix(n, 40);

    std::vector<int> ks = {1, 10, 100};

    for (int k : ks) {
        std::vector<Vector> rhs(k);
        for (int i = 0; i < k; i++)
            rhs[i] = generate_random_vector(n, 200 + i);

        Timer t;

        //гаусс: каждый раз полный алгоритм заново
        t.begin();
        for (int i = 0; i < k; i++)
            gauss_pivot(A, rhs[i]);
        double t_gauss = t.elapsed_ms();

        //LU: разложение один раз, подстановки k раз
        t.begin();
        LU lu = lu_decompose(A);
        for (int i = 0; i < k; i++)
            lu_solve(lu, rhs[i]);
        double t_lu = t.elapsed_ms();

        std::cout << k << "\t" << t_gauss << "\t\t\t" << t_lu << "\n";
    }
}

//эксперимент 3: матрица гильберта — плохо обусловленные системы
void experiment_3() {
    std::cout << "\n=== Experiment 3: Hilbert matrix (ill-conditioned) ===\n";
    std::cout << "n\tmethod\t\terror\t\t\tresidual\n";
    std::cout << "------------------------------------------------------------\n";

    std::vector<int> sizes = {5, 10, 15};

    for (int n : sizes) {
        Matrix H = generate_hilbert(n);

        //точное решение x = (1, 1, ..., 1)
        Vector x_exact(n, 1.0);

        //правая часть b = H * x_exact
        Vector b(n, 0.0);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                b[i] += H[i][j] * x_exact[j];

        //гаусс без выбора
        try {
            Vector x = gauss_no_pivot(H, b);
            std::cout << n << "\tno pivot\t" << relative_error(x, x_exact) << "\t\t" << residual(H, x, b) << "\n";
        } catch (...) {
            std::cout << n << "\tno pivot\terror\n";
        }

        //гаусс с выбором
        try {
            Vector x = gauss_pivot(H, b);
            std::cout << n << "\tpivot\t\t" << relative_error(x, x_exact) << "\t\t" << residual(H, x, b) << "\n";
        } catch (...) {
            std::cout << n << "\tpivot\t\terror\n";
        }

        //LU
        try {
            LU lu = lu_decompose(H);
            Vector x = lu_solve(lu, b);
            std::cout << n << "\tLU\t\t" << relative_error(x, x_exact) << "\t\t" << residual(H, x, b) << "\n";
        } catch (...) {
            std::cout << n << "\tLU\t\tsingular\n";
        }
    }
}

int main() {
    experiment_1();
    experiment_2();
    experiment_3();
    return 0;
}