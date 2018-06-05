#include "SOM.hpp"

SomIndex::SomIndex(int inX = 0, int inY = 0)
{
	x = inX;
	y = inY;
}

SomIndex::~SomIndex()
{
	;
}
/*
int SomIndex::getSomIndex(Som map)
{
	return map.getWidth()*y + x;
}
*/

unsigned int SomIndex::getX() const
{
	return x;
}

unsigned int SomIndex::getY() const
{
	return y;
}

void SomIndex::setX(int ix)
{
	x = ix;
}

void SomIndex::setY(int iy)
{
	y = iy;
}