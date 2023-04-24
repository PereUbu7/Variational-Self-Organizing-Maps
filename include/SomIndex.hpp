#pragma once

#include <stddef.h>

class Som;

class SomIndex
{
protected:
	size_t x, y;

public:
	SomIndex(size_t x, size_t y) noexcept;
	SomIndex(const Som &map, size_t index) noexcept;
	~SomIndex() = default;
	size_t getSomIndex(const Som &som);
	size_t getX() const noexcept;
	size_t getY() const noexcept;
	void setX(size_t index) noexcept;
	void setY(size_t index) noexcept;

	bool operator==(const SomIndex &other)
	{
		return x == other.x &&
			y == other.y;
	}
};