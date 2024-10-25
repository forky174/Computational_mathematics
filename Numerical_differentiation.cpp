#include <iostream>
#include <array>
#include <cmath>

template<typename RealType, unsigned long int N>
struct DerivativeCoef {
    RealType centralCoef;
    std::array<RealType, N> otherCoefs;
};

template<typename RealType, unsigned long int N>
DerivativeCoef<RealType, N> calcDerivativeCoef(const std::array<RealType, N>& points) noexcept {
    //create matrix for linear system
    std::array<std::array<RealType, N+2>, N+1> matrix;
    for (unsigned long int i = 0; i < N+1; ++i)
        std::fill(std::begin(matrix[i]), std::begin(matrix[i]) + N + 2, 0);

    for (unsigned long int i = 0; i < N+1; ++i) {
        for (unsigned long int j = 0; j < N+1; ++j) {
            if (i == 0)
                matrix[i][j] = 1;
            else {
                if (j == 0)
                    matrix[i][j] = 0;
                else
                    matrix[i][j] = std::pow((points[j-1]),i) / std::tgamma(i+1);
            } 
        }
    }
    matrix[1][N+1] = 1; // задаем искомый порядок производной

    // std::cout << "Матрица:" << std::endl;
    // for (unsigned long int i = 0; i < N+1; ++i) {
    //     for (unsigned long int j = 0; j < N+2; ++j)
    //         std::cout << matrix[i][j] << " ";
    //     std::cout << std::endl;
    // }

    // прямой ход
    RealType tmp;
    std::array<RealType, N+1> Coefs;

    for (unsigned long int k = 1; k < N+1; k++) {
        for (unsigned long int j = k; j < N+1; j++) {
            tmp = matrix[j][k-1] / matrix[k-1][k-1];
            for (unsigned long int i = 0; i < N+2; i++)
                matrix[j][i] = matrix[j][i] - tmp*matrix[k-1][i];
        }
    }

    // обратный ход
    for (unsigned long int i = N; (i >= 0) && (i < N+1); i--) {
        Coefs[i] = matrix[i][N+1] / matrix[i][i];
        for (unsigned long int c = N; (c > i) && (c < N+1); c--) {
            Coefs[i] = Coefs[i] -  matrix[i][c] * Coefs[c] / matrix[i][i];
        }
    }

    DerivativeCoef<RealType, N> result;
    result.centralCoef = Coefs[0];
    for (unsigned long int i = 0; i < N; ++i)
        result.otherCoefs[i] = Coefs[i+1];
    return result;
};

float fun (float x) {
    return std::exp(x);
}

int main () {
    const long unsigned int N = 2;
    std::array<long double, N> points = {-1, 1};

    // const long unsigned int N = 3;
    // std::array<float, N> points = {-1, 1, 2};

    // const long unsigned int N = 4;
    // std::array<float, N> points = {-2, -1, 1, 2};

    // const long unsigned int N = 5;
    // std::array<float, N> points = {-2, -1, 1, 2, 3};

    std::array<long double, 20> hs = {4, 3, 2, 1, 0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15};
    long double x0 = 1.;
    std::array<long double, 20> defs;

    for (long unsigned int i = 0; i < hs.size(); ++i) {
        std::array<long double, N> points_h;
        for (unsigned long int j = 0; j < N; ++j)
            points_h[j] = points[j] * hs[i];
        DerivativeCoef<long double, N> coefs = calcDerivativeCoef(points_h);
        defs[i] = coefs.centralCoef * fun(x0);
        for (long unsigned int n = 0; n < N; ++n) {
            defs[i] += coefs.otherCoefs[n] * fun(x0 + points_h[n]);
        }
        // std::cout << "Значение производной при h=" << hs[i] << " :: " << defs[i] << std::endl;
        std::cout << "Значение ошибки      при h=" << hs[i] << " :: " << std::abs(fun(x0)-defs[i]) << std::endl;
    }
}