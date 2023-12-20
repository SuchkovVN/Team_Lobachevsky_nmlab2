#include "nmlib.hpp"
#include "nmlinalg.hpp"
#include <iostream>
#include "diffschemes.hpp"

float test_rhs(float x, float v, float y) {
    return -(4. / 2.) * v;
}

float test_rhs2(float x, float v, float y) {
    return v;
}

int main() {

    double coefs[6][3] = {
        {0.l, -0.7367, 1.l},
        {1.l, -0.65573, -0.06059},
        {-0.06059, -0.49614, 1.l},
        {1.l, -0.78194, -0.1881},
        {-0.17941, -0.78194, 1.l},
        {1.l, -0.64833, 0.l}
    };

    double rhs[6] = {0.520632, 0.0278, 0.0196, 0.03582, 0.34792, 0.38185};

    double vars[6] = {0.l, 0.l, 0.l, 0.l, 0.l, 0.l};

    TridiagonalAlg(coefs, vars, rhs, 6);

    for(size_t i = 0; i < 6; i++) {
        std::cout << vars[i] << '\n';
    }

    return 0;
}