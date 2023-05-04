#pragma once

#include <string>

struct ColumnSpec
{
    ColumnSpec(const std::string &name, const float &weight, const int &isBinary) :
        name{name},
        weight{weight},
        isBinary{isBinary} {}
    ColumnSpec(const std::string &&name, const float &&weight, const int &&isBinary) = delete;
    ColumnSpec(std::string &&name, float &&weight, int &&isBinary) = delete;
    
    const std::string &name;
    const float &weight;
    const int &isBinary;
};