#pragma once

#include <cstddef>
template <class T, class U>
void TridiagonalAlg(T coeffs, U x, U rhs, const size_t& n) {
    /* coeffs is an object, which contains a TMA coefficients 
       and can that coefficients can be accessed in form coeffs[i][j] where j =0,1,2 
       x is an object, which represent a vector of linear system variables which can 
       be accessed in form x[i] (x[0] = x0, x[1] = x1...)
       n is number of variables in linear system 
       rhs is an object, which represent a vector of linear system right-hand-sides which can 
       be accessed in form rhs[i] (rhs[0] = F1 (?)) */

    double a[2] = {0.l, 0.l};
    double b[2] = {0.l, 0.l};


    for (size_t i = 0; i < n; i++) {
        a[1] = -coeffs[i][2] / (coeffs[i][0] * a[0] + coeffs[i][1]);
        b[1] = (rhs[i] - coeffs[i][0] * b[0]) / (coeffs[i][0] * a[0] + coeffs[i][1]);

        x[i + 1] = (x[i] - b[1]) / a[1]; 
        a[0] = a[1];
        b[0] = b[1];
    }
}