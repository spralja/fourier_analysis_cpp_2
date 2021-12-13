//
// Created by Robert on 13/12/2021.
//

#ifndef FOURIER_ANALYSIS_CPP_2_COEFFICIENTCOLLECTION_H
#define FOURIER_ANALYSIS_CPP_2_COEFFICIENTCOLLECTION_H

#include <unordered_map>
#include <memory>

class FourierAnalysis;

class CoefficientCollection {
private:
    const FourierAnalysis* parent;
    mutable std::unordered_map<int, std::unordered_map<int, std::unordered_map<int, std::unique_ptr<double>>>>
        coefficients;

public:
    explicit CoefficientCollection(const FourierAnalysis* parent);

    const double& get(const int& k, const int& n, const int& m) const;
};


#endif //FOURIER_ANALYSIS_CPP_2_COEFFICIENTCOLLECTION_H
