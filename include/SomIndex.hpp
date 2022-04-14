#pragma once

#include <stddef.h>

class SomIndex
{
	protected:
		size_t x,y;
	
	public:
		SomIndex(size_t, size_t) noexcept;
		~SomIndex() = default;
		//int getSomIndex(Som);
		size_t getX() const noexcept;
		size_t getY() const noexcept;
		void setX(size_t) noexcept;
		void setY(size_t) noexcept;
};