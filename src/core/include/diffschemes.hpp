#pragma once

#include "Table.hpp"
#include "logger.hpp"
#include "nmlib.hpp"

#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <vector>

struct Uniform1DNet {
    /* Uniform1DNet is uniform net with nodes: 
    xi = first + i * step for i=0,n and step = (last - first / n) */

    double first;
    double last;
    size_t n;

    double step;

    Uniform1DNet(const double& _a, const double& _b, const double& _n) {
        first = _a;
        last = _b;
        n = _n;

        step = (last - first) / n;
    }
};

class balanceMethod {
    using func = std::function<double(double, double)>;
    std::shared_ptr<Uniform1DNet> _net;
    std::unique_ptr<Table> _result;

    /* coeff functions for diff scheme */
    func _ck;
    func _cq;
    func _cf;

public:
    balanceMethod();
    virtual ~balanceMethod() {
        _net.reset();
    }

    balanceMethod(Uniform1DNet* net, Table* tbl, func&& ck, func&& cq, func&& cf) : _net(net), _result(tbl), _ck(ck), _cq(cq), _cf(cf) {}

    virtual void eval();

    virtual Table* getTable() {
        return _result.release();
    }

    virtual void setNet(Uniform1DNet* net) {
        _net.reset(net);
    }
};

class NMbalance : public balanceMethod {
using func = std::function<double(double, double)>;
public:
    NMbalance();
    ~NMbalance() override;

    NMbalance(Uniform1DNet* net, Table* tbl, func&& ck, func&& cq, func&& cf) : balanceMethod(net, tbl, std::move(ck), std::move(cq), std::move(cf)) {}

    void eval() override;
private:
    std::vector<double> matrix; // matrix of linear equations system, with C-style layout
};