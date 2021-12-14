//
// Created by Robert on 13/12/2021.
//

#ifndef FOURIER_ANALYSIS_CPP_2_TRIGONOMETRICFUNCTIONS_H
#define FOURIER_ANALYSIS_CPP_2_TRIGONOMETRICFUNCTIONS_H

#include <vector>

class FourierAnalysis;

class TrigonometricFunctions {
private:
    const FourierAnalysis* parent;

    std::vector<double> sinPhiValues, cosPhiValues, tanPhiValues, sinThetaValues, cosThetaValues, tanThetaValues;
public:
    explicit TrigonometricFunctions(const FourierAnalysis* parent);

    [[nodiscard]] const double& sinPhi(const int& index) const;
    [[nodiscard]] const double& cosPhi(const int& index) const;
    [[nodiscard]] const double& tanPhi(const int& index) const;
    [[nodiscard]] const double& sinTheta(const int& index) const;
    [[nodiscard]] const double& cosTheta(const int& index) const;
    [[nodiscard]] const double& tanTheta(const int& index) const;
};


#endif //FOURIER_ANALYSIS_CPP_2_TRIGONOMETRICFUNCTIONS_H
