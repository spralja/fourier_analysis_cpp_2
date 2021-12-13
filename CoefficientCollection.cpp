//
// Created by Robert on 13/12/2021.
//

#include "CoefficientCollection.h"

#include <memory>
#include "FourierAnalysis.h"

CoefficientCollection::CoefficientCollection(const FourierAnalysis *parent): parent(parent) {}

const double &CoefficientCollection::get(const int &k, const int &n, const int &m) const {
    if(coefficients[k][n][m] == nullptr)
        coefficients[k][n][m] = std::make_unique<double>(parent->C(k, n, m));

    return *coefficients[k][n][m];
}
