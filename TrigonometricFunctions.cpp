//
// Created by Robert on 13/12/2021.
//

#include "TrigonometricFunctions.h"
#include "FourierAnalysis.h"

TrigonometricFunctions::TrigonometricFunctions(const FourierAnalysis *parent): parent(parent) {
    sinPhiValues.reserve(parent->n_phi);
    cosPhiValues.reserve(parent->n_phi);
    tanPhiValues.reserve(parent->n_phi);
    sinThetaValues.reserve(parent->n_theta);
    cosThetaValues.reserve(parent->n_theta);
    tanThetaValues.reserve(parent->n_theta);
    for(int phi_index = 0; phi_index < parent->n_phi; ++phi_index) {
        sinPhiValues[phi_index] = std::sin(parent->phi(phi_index));
        cosPhiValues[phi_index] = std::cos(parent->phi(phi_index));
        tanPhiValues[phi_index] = std::tan(parent->phi(phi_index));
    }

    for(int theta_index = 0; theta_index < parent->n_theta; ++theta_index) {
        sinThetaValues[theta_index] = std::sin(parent->theta(theta_index));
        cosThetaValues[theta_index] = std::cos(parent->theta(theta_index));
        tanThetaValues[theta_index] = std::tan(parent->theta(theta_index));
    }
}

const double &TrigonometricFunctions::sinPhi(const int &index) const {
    return sinPhiValues[index];
}

const double &TrigonometricFunctions::cosPhi(const int &index) const {
    return cosPhiValues[index];
}

const double &TrigonometricFunctions::tanPhi(const int &index) const {
    return tanPhiValues[index];
}

const double &TrigonometricFunctions::sinTheta(const int &index) const {
    return sinThetaValues[index];
}

const double &TrigonometricFunctions::cosTheta(const int &index) const {
    return cosThetaValues[index];
}

const double &TrigonometricFunctions::tanTheta(const int &index) const {
    return tanThetaValues[index];
}

