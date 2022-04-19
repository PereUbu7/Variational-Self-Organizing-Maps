#pragma once

#include "SomIndex.hpp"
#include <vector>

class UMatrix
{
private:
    std::vector<double> data;
    size_t width, height;
public:
    UMatrix(std::vector<double> Data, size_t Width, size_t Height) :
        data{Data},
        width{Width},
        height{Height}
    {}
    double getValueAtIndex(size_t Width, size_t Height) const;
    double getValueAtIndex(SomIndex Index) const;
    const std::vector<double>& getData() const noexcept;
    size_t getWidth() const noexcept;
    size_t getHeight() const noexcept;
    double getMinValue() const noexcept;
    double getMaxValue() const noexcept;
};