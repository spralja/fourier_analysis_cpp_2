//
// Created by Robert on 13/12/2021.
//

#include "TrigonometricFunctions.h"
#include "FourierAnalysis.h"

TrigonometricFunctions::TrigonometricFunctions(const FourierAnalysis *parent): parent(parent) {
    const int& n_phi = parent->n_phi;
    const int& n_theta = parent->n_theta;
    for(int i = 0; i < n_phi; ++i) {
        const double phi = parent->phi(i);
        sinValues[phi] = std::sin(phi);
        cosValues[phi] = std::cos(phi);
        tanValues[phi] = std::tan(phi);
    }

    for(int i = 0; i < n_theta; ++i) {
        const double theta = parent->theta(i);
        sinValues[theta] = std::sin(theta);
        cosValues[theta] = std::cos(theta);
        tanValues[theta] = std::tan(theta);
    }
}

const double& TrigonometricFunctions::sin(const double &beta) const {
    return sinValues[beta];
}

const double &TrigonometricFunctions::cos(const double &beta) const {
    return cosValues[beta];
}

const double &TrigonometricFunctions::tan(const double &beta) const {
    return tanValues[beta];
}




