#include "diffschemes.hpp"

#include "Table.hpp"
#include "logger.hpp"
#include "nmlib.hpp"
#include "nmlinalg.hpp"
#include <cmath>
#include <cstddef>
#include <math.h>
#include <vector>

NMbalance::~NMbalance() {}

void NMbalance::eval() {
    const double x = _net->first;  // x coord of node in net
    double step = _net->step;
    size_t n = _net->n;
    double stepsq = step * step;

    LOG_DEBUG_CLI("net: ", _net->n, _net->step, _net->first, _net->last);
    LOG_DEBUG_CLI("stepsq:", stepsq);
    // evaluating the coefficients for linear system (tridiagonal matrix algorithm(TMA) aka "progonka")
    auto cA = [&](const size_t& i) -> double {
        return _ca(x + (i - 1) * step, x + i * step, step) / stepsq;
    };
    auto cB = [&](const size_t& i) -> double {
        return -((_ca(x + (i - 1) * step, x + i * step, step) + _ca(x + i * step, x + (i + 1) * step, step)) / stepsq +
                 _cd(x + (i - 0.5) * step, x + (i + 0.5) * step, step));
    };
    auto cC = [&](const size_t& i) -> double {
        return _ca(x + i * step, x + (i + 1) * step, step) / stepsq;
    };

    auto cF = [&](const size_t& i) -> double {
        return -_cphi(x + (i - 0.5) * step, x + (i + 0.5) * step, step);
    };

    LOG_DEBUG_CLI("Evalutating matrix coeffs...");
    matrix[0][0] = 0.l;
    matrix[0][1] = 1.l;
    matrix[0][2] = 0.l;
    rhs[0] = _mu1;

    for (size_t i = 1; i < n; i++) {
        matrix[i][0] = cA(i);
        matrix[i][1] = cB(i);
        matrix[i][2] = cC(i);
        rhs[i] = cF(i);
    }

    matrix[n][0] = 0.l;
    matrix[n][1] = 1.l;
    matrix[n][2] = 0.l;
    rhs[n] = _mu2;
    // coefficient evaluated and wrote in matrix and rhs vars
    // now we need to solve linear system w/ TMA

    LOG_DEBUG_CLI("Applying TMA algorithm...");
    TridiagonalAlg(matrix, vars, rhs, n + 1);

    for (size_t i = 0; i < vars.size(); i++) {
        double point_x = x + step * i;
        double sol_value = _sol(point_x);
        _result->append({ point_x, vars[i], sol_value, std::abs(vars[i] - sol_value), 0.l });
    }

    // half-step calc

    step = step / 2;
    n = n * 2;
    stepsq = stepsq / 4;
    LOG_DEBUG_CLI("half-step", step, stepsq, n);
    std::vector<std::vector<double>> matrix2{};
    std::vector<double> vars2(n);
    std::vector<double> rhs2(n);
    matrix2.resize(n + 1);
    for (auto& m : matrix2) {
        m.resize(3);
    }

    matrix2[0][0] = 0.l;
    matrix2[0][1] = 1.l;
    matrix2[0][2] = 0.l;
    rhs2[0] = _mu1;

    for (size_t i = 1; i < n; i++) {
        matrix2[i][0] = cA(i);
        matrix2[i][1] = cB(i);
        matrix2[i][2] = cC(i);
        rhs2[i] = cF(i);
    }

    matrix2[n][0] = 0.l;
    matrix2[n][1] = 1.l;
    matrix2[n][2] = 0.l;
    rhs2[n] = _mu2;
    // coefficient evaluated and wrote in matrix and rhs vars
    // now we need to solve linear system w/ TMA

    LOG_DEBUG_CLI("Applying TMA algorithm...");
    TridiagonalAlg(matrix2, vars2, rhs2, n + 1);

    for (size_t i = 0; i < vars.size(); i++) {
        (*_result)[i].v_2 = vars2[2 * i];
        (*_result)[i].vdiff = std::abs((*_result)[i].v - (*_result)[i].v_2);
    }
    // n = n / 2;
    // step = step * 2;
}