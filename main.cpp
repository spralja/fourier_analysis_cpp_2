#include <iostream>
#include "FourierAnalysis.h"

int main() {
    auto fa = FourierAnalysis(1000, 1000, 1);
    std::cout << fa.C(1, 1, 1);
    std::cout << fa.fourierSum(1, 1, 1);
    return 0;
}
