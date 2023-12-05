#pragma once

#include <vector>
struct TableRow {

};

class Table {
using rows = std::vector<TableRow>;
    rows _data;
public:
    Table();
    virtual ~Table();

    Table(rows&& data) : _data(data) {}
};