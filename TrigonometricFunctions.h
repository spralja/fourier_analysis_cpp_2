//
// Created by Robert on 13/12/2021.
//

#ifndef FOURIER_ANALYSIS_CPP_2_TRIGONOMETRICFUNCTIONS_H
#define FOURIER_ANALYSIS_CPP_2_TRIGONOMETRICFUNCTIONS_H

#include <unordered_map>

class FourierAnalysis;

class TrigonometricFunctions {
private:
    const FourierAnalysis* parent;

    mutable std::unordered_map<double, double> sinValues;
    mutable std::unordered_map<double, double> cosValues;
    mutable std::unordered_map<double, double> tanValues;
public:
    explicit TrigonometricFunctions(const FourierAnalysis* parent);

    const double& sin(const double& beta) const;
    const double& cos(const double& beta) const;
    const double& tan(const double& beta) const;
};


#endif //FOURIER_ANALYSIS_CPP_2_TRIGONOMETRICFUNCTIONS_H
