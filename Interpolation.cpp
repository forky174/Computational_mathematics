#include <array>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <typeinfo>
#define _USE_MATH_DEFINES

template<typename xType, typename yType, unsigned int N>
class NewtonInterpolation {
    private:
    std::array<xType, N> xs;
    std::array<yType, N> div_difs;
    public:
    void NewtonInterpolator(const std::array<xType, N> &points, const std::array<yType, N>& values) noexcept {
        for (unsigned long int n = 0; n < N; ++n) {
            xs[n] = points[n];
            div_difs[n] = 0.;
            for (unsigned long int j = 0; j <= n; ++j) {
                long double denominator = 1.;
                for (unsigned long int i = 0; i <= n ; ++i)
                    if (i != j)
                        denominator *= (points[j] - points[i]);
                div_difs[n] += values[j] / denominator;
            }
        }
    };
    yType interpolate(const xType& x) const noexcept {
        yType result = div_difs[N-1];
        for (unsigned long int k = N-2; (k >= 0) && (k < N-1); --k)
            result = div_difs[k] + (x - xs[k]) * result;
        return result;
    };
};

long double fun(long double x) {
    return std::exp(x);
};

int main () {
    const long unsigned int N = 5;
    std::array<std::pair<long double, long double>, 6> line_segments {std::make_pair(0., 2.),    std::make_pair(0., 1.),
                                                                      std::make_pair(0., 1./2.), std::make_pair(0., 1./4.),
                                                                      std::make_pair(0., 1./8.), std::make_pair(0., 1./16.)};


    // построение равномерного распределения узлов для каждого отрезка
    std::array<std::array<long double, N>, 6> Uniform_distribution_nodes;
    for (unsigned long int i = 0; i < line_segments.size(); ++i) {
        long double step = (line_segments[i].second - line_segments[i].first) / (N-1);
        // std::cout << "Отрезок " <<  line_segments[i].second - line_segments[i].first << std::endl;
        for (unsigned long j = 0; j < N; ++j) {
            Uniform_distribution_nodes[i][j] = line_segments[i].first + step * j;
            // std::cout << Uniform_distribution_nodes[i][j] << " ";
        }
        // std::cout << std::endl;
    }

    std::array<std::array<long double, N>, 6> Uniform_distribution_values;
    for (unsigned long int i = 0; i < line_segments.size(); ++i)
        for (unsigned long j = 0; j < N; ++j) {
            Uniform_distribution_values[i][j] = fun(Uniform_distribution_nodes[i][j]);
            // std::cout << Uniform_distribution_values[i][j] << " ";
        }
        // std::cout << std::endl;

    // построение Чебышевского распределения узлов для каждого отрезка
    std::array<std::array<long double, N>, 6> Chebyshev_nodes;
    for (unsigned long int i = 0; i < line_segments.size(); ++i) {
        long double a = line_segments[i].first;
        long double b = line_segments[i].second;
        // std::cout << "Отрезок " <<  b - a << std::endl;
        for (unsigned long j = N; j >= 1; j--) {
            Chebyshev_nodes[i][j-1] = (a+b)/2 + (b-a)/2 * std::cos(M_PI * (2*j - 1) / (2*N));
            // std::cout << Chebyshev_nodes[i][j-1] << " ";
        }
        // std::cout << std::endl;
    }

    std::array<std::array<long double, N>, 6> Chebyshev_values;
    for (unsigned long int i = 0; i < line_segments.size(); ++i)
        for (unsigned long j = 0; j < N; ++j)
            Chebyshev_values[i][j] = fun(Chebyshev_nodes[i][j]);

    // рассчет интерполянта Ньютона и получение ошибок
    for (unsigned long int i = 0; i < line_segments.size(); ++i) {
        long double length = line_segments[i].second - line_segments[i].first;
        long double x = line_segments[i].first;
        long double step = length / 999;
        long double Uniform_distribution_error = 0;
        long double Chebyshev_error = 0;
        long double Uniform_distribution_inter_value;
        long double Chebyshev_inter_value;
        
        NewtonInterpolation<long double, long double, N> Uniform_distribution;
        NewtonInterpolation<long double, long double, N> Chebyshev;

        // std::reverse(Chebyshev_nodes[i].begin(), Chebyshev_nodes[i].end());
        // std::reverse(Chebyshev_values[i].begin(), Chebyshev_values[i].end());

        Uniform_distribution.NewtonInterpolator(Uniform_distribution_nodes[i], Uniform_distribution_values[i]);
        Chebyshev.NewtonInterpolator(Chebyshev_nodes[i], Chebyshev_values[i]);

        for (unsigned long int p = 0; p < 1000; ++p) {
            // std::cout << x << std::endl;
            Uniform_distribution_inter_value = Uniform_distribution.interpolate(x);
            Chebyshev_inter_value = Chebyshev.interpolate(x);
            // std::cout << "Равномерное распределение на отрезке " << length << " :: " << Uniform_distribution_inter_value << std::endl;
            // std::cout << "Чебышевское распределение на отрезке " << length << " :: " << Chebyshev_inter_value << std::endl;
            // std::cout << "Значения функции на отрезке          " << length << " :: " << fun(x) << std::endl;
            if (std::abs(Uniform_distribution_inter_value - fun(x)) > Uniform_distribution_error)
                Uniform_distribution_error  = std::abs(Uniform_distribution_inter_value - fun(x));
            if (std::abs(Chebyshev_inter_value - fun(x)) > Chebyshev_error)
                Chebyshev_error  = std::abs(Chebyshev_inter_value - fun(x));
            x += step;
        }
        std::cout << "Ошибка при равномерном распределении узлов при длине отрезка " << length << " :: " << Uniform_distribution_error << std::endl;
        std::cout << "Ошибка при Чебышевском распределении узлов при длине отрезка " << length << " :: " << Chebyshev_error << std::endl;
    }
}