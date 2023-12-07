#include "diffschemes.hpp"

#include "nmlinalg.hpp"
#include "logger.hpp"
#include <cmath>
#include <cstddef>

void NMbalance::eval() {
    const double& x = _net->first;  // x coord of node in net
    const double& step = _net->step;
    const double& n = _net->n;
    const double stepsq = step * step;

    // evaluating the coefficients for linear system (tridiagonal matrix algorithm(TMA) aka "progonka")
    auto cA = [&](const size_t& i) -> double {
        return _ca(x + (i - 1) * step, x + i * step) / stepsq;
    };
    auto cB = [&](const size_t& i) -> double {
        return -(_ca(x + (i - 1) * step, x + i * step) / stepsq + _ca(x + (i)*step, x + (i + 1) * step) / stepsq +
                 _cd(x + (i - 0.5) * step, x + (i + 0.5) * step));
    };
    auto cC = [&](const size_t& i) -> double {
        return _ca(x + i * step, x + (i + 1) * step) / stepsq;
    };

    auto cF = [&](const size_t& i) -> double {
        return -_cphi(x + (i - 0.5) * step, x + (i + 0.5) * step);
    };

    matrix[0][0] = 0.l;
    matrix[0][1] = cB(1);
    matrix[0][2] = cC(1);
    rhs[0] = cF(1) - cA(1) * _mu1;

    for (size_t i = 1; i < n; i++) {
        matrix[i][0] = cA(i + 1);
        matrix[i][1] = cB(i + 1);
        matrix[i][2] = cC(i + 1);
        rhs[i] = cF(i + 1);
    }

    matrix[n][0] = cA(n + 1);
    matrix[n][1] = cB(n + 1);
    matrix[n][2] = 0.l;
    rhs[n] = cF(n + 1) - cC(n + 1) * _mu2;
    // coefficient evaluated and wrote in matrix and rhs vars
    // now we need to solve linear system w/ TMA

    TridiagonalAlg(matrix, vars, rhs, n + 1);

    NM_ASSERT((std::abs(vars[0] -_mu1) < 1e-9) && (std::abs(vars[n] - _mu2) < 1e-9), "error while eval: incorrect value");

    // actually we need to put solution in some resultTable like object, but for some time it will be more... reasonable
}