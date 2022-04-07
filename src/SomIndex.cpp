#include "SomIndex.hpp"

SomIndex::SomIndex(size_t inX = 0, size_t inY = 0) noexcept :
	x{inX},
	y{inY}
{}

/*
int SomIndex::getSomIndex(Som map)
{
	return map.getWidth()*y + x;
}
*/

size_t SomIndex::getX() const noexcept
{
	return x;
}

size_t SomIndex::getY() const noexcept
{
	return y;
}

void SomIndex::setX(size_t ix) noexcept
{
	x = ix;
}

void SomIndex::setY(size_t iy) noexcept
{
	y = iy;
}