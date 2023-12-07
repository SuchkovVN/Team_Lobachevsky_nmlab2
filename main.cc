#include "nmlib.hpp"
#include "nmlinalg.hpp"
#include <iostream>

float test_rhs(float x, float v, float y) {
    return -(4. / 2.) * v;
}

float test_rhs2(float x, float v, float y) {
    return v;
}

int main() {
    double matrix[3][3] = {{0., 2., 2.,}, {1.,1.5,2.}, {1., 5., 0.,}};
    double rhs[3] = {0.5, 0.8, 0.9};
    double v[3] = {1., 0., 1.,};

    TridiagonalAlg(matrix, v, rhs, 3);
    std::cout << v[0] << '\n' << v[1] << '\n' << v[2] << '\n';
    return 0;
}