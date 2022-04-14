#pragma once

#include "SomIndex.hpp"
#include <vector>

template <typename Precision>
class UMatrix
{
private:
    std::vector<Precision> data;
    size_t width, height;
public:
    UMatrix(std::vector<Precision> Data, size_t Width, size_t Height) :
        data{Data},
        width{Width},
        height{Height}
    {}
    Precision getValueAtIndex(size_t Width, size_t Height) const;
    Precision getValueAtIndex(SomIndex Index) const;
    const std::vector<Precision>& getData() const noexcept;
    size_t getWidth() const noexcept;
    size_t getHeight() const noexcept;
};