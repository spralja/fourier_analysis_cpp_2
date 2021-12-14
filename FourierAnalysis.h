//
// Created by Robert on 13/12/2021.
//

#ifndef FOURIER_ANALYSIS_CPP_2_FOURIERANALYSIS_H
#define FOURIER_ANALYSIS_CPP_2_FOURIERANALYSIS_H


#include <cmath>
#include <unordered_map>

#include "TrigonometricFunctions.h"
#include "CoefficientCollection.h"

#include <vector>

class FourierAnalysis {
public:
    const double a_phi, b_phi, a_theta, b_theta;
    const int n_phi, n_theta, n_sigma;
    const double d_phi, d_theta;
    const double mu;
private:
    const TrigonometricFunctions trigs;
    const CoefficientCollection coefficients;
public:
    FourierAnalysis(const int& n_phi, const int& n_theta, const int& n_sigma, const double& a_phi = -M_PI / 2,
                    const double& b_phi = M_PI / 2, const double& a_theta = 0, const double& b_theta = 2 * M_PI,
                    const double& mu = -1 / (std::sqrt(2 * M_PI) * std::sqrt(2 * M_PI) * std::sqrt(2 * M_PI)));

    double phi(const int& index) const;

    double theta(const int& index) const;

    double C(const int& k, const int& n, const int& m) const;

    double getC(const int& k, const int& n, const int& m) const;

    double fourierSum(const double& x, const double& y, const double& z) const;

    friend class TrigonometricFunctions;
};


#endif //FOURIER_ANALYSIS_CPP_2_FOURIERANALYSIS_H
