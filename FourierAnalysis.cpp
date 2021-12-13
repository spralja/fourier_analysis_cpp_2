//
// Created by Robert on 13/12/2021.
//

#include "FourierAnalysis.h"

FourierAnalysis::FourierAnalysis(const int &n_phi, const int &n_theta, const int& n_sigma, const double &a_phi, const double &b_phi,
                                 const double &a_theta, const double &b_theta)
                                 : n_phi(n_phi), n_theta(n_theta), n_sigma(n_sigma), a_phi(a_phi), b_phi(b_phi),
                                 a_theta(a_theta), b_theta(b_theta), d_phi((b_phi - a_phi) / n_phi),
                                 d_theta((b_theta - a_theta) / n_theta), trigs(this), coefficients(this) {}

double FourierAnalysis::phi(const int& index) const {
    return a_phi + d_phi * (index + 0.5);
}

double FourierAnalysis::theta(const int& index) const {
    return a_theta + d_theta * (index + 0.5);
}

double FourierAnalysis::C(const int& k, const int& n, const int& m) const {
    double sum = 0.0;
    for(int i = 0; i < n_phi; ++i)
        for(int j = 0; j < n_theta; ++j)
            sum += 1 / (
                    m * m * trigs.sin(phi(i)) * trigs.tan(phi(i)) +
                    2 * m * trigs.sin(phi(i)) * (k * trigs.cos(theta(j)) + n * trigs.sin(theta(j))) +
                    trigs.cos(phi(i)) *
                    (k * trigs.cos(theta(j)) + n * trigs.sin(theta(j))) *
                            (k * trigs.cos(theta(j)) + n * trigs.sin(theta(j)))
                    );

    sum *= -1 / (std::sqrt(2 * M_PI) * std::sqrt(2 * M_PI) * std::sqrt(2 * M_PI)) * d_phi * d_theta;
    return sum;
}

double FourierAnalysis::fourierSum(const double &x, const double &y, const double &z) const {
    double sum = 0.0;
    for(int k = -n_sigma; k <= n_sigma; ++k)
        for(int n = -n_sigma; n <= n_sigma; ++n)
            for(int m = 1; m <= n_sigma; ++m)
                sum += coefficients.get(k, n, m) * std::cos(k * x + n * y + m * z);

    for(int k = -n_sigma; k <= n_sigma; ++k)
        for(int n = 1; n <= n_sigma; ++n)
            sum += coefficients.get(k, n, 0) * std::cos(k * x + n * y);

    for(int k = 1; k <= n_sigma; ++k)
        sum += coefficients.get(k, 0, 0) * std::cos(k * x);

    return 2 * sum;
}




