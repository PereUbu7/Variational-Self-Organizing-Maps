#include "UMatrix.hpp"

template<typename Precision>
Precision UMatrix<Precision>::getValueAtIndex(size_t inWidth, size_t inHeight) const
{
    auto index = width*inHeight + inWidth;
    assert(index <= data.size());
    return data[index];
}

template<typename Precision>
Precision UMatrix<Precision>::getValueAtIndex(SomIndex index) const
{
    auto realIndex = width*index.getY() + index.getX();
    assert(realIndex <= data.size());
    return data[realIndex];
}

template<typename Precision>
const std::vector<Precision>& UMatrix<Precision>::getData() const noexcept
{
    return data;
}

template<typename Precision>
size_t UMatrix<Precision>::getWidth() const noexcept
{
    return width;
}

template<typename Precision>
size_t UMatrix<Precision>::getHeight() const noexcept
{
    return height;
}