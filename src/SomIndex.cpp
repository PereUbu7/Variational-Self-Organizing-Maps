#include "SomIndex.hpp"

class Som
{
	public: 
		size_t getWidth() const noexcept;
		size_t getHeight() const noexcept;
};

SomIndex::SomIndex(size_t inX = 0, size_t inY = 0) noexcept :
	x{inX},
	y{inY}
{}

SomIndex::SomIndex(const Som &map, size_t index) noexcept :
	x{index % map.getWidth()},
	y{(index - index % map.getWidth())/map.getHeight()}
	{}


size_t SomIndex::getSomIndex(const Som &map)
{
	return map.getWidth()*y + x;
}


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