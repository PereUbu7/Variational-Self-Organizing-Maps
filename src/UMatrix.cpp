#include "UMatrix.hpp"

#include <algorithm>
#include <cassert>

double UMatrix::getValueAtIndex(size_t inWidth, size_t inHeight) const
{
    auto index = width*inHeight + inWidth;
    assert(index <= data.size());
    return data[index];
}

double UMatrix::getValueAtIndex(SomIndex index) const
{
    auto realIndex = width*index.getY() + index.getX();
    assert(realIndex <= data.size());
    return data[realIndex];
}

const std::vector<double>& UMatrix::getData() const noexcept
{
    return data;
}

size_t UMatrix::getWidth() const noexcept
{
    return width;
}

size_t UMatrix::getHeight() const noexcept
{
    return height;
}

double UMatrix::getMinValue() const noexcept
{
    return *std::min_element(data.begin(), data.end());
}

double UMatrix::getMaxValue() const noexcept
{
    return *std::max_element(data.begin(), data.end());
}