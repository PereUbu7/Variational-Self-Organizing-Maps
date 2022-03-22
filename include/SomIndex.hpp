#include <stddef.h>

class SomIndex
{
	protected:
		size_t x,y;
	
	public:
		SomIndex(size_t, size_t);
		~SomIndex() = default;
		//int getSomIndex(Som);
		size_t getX() const;
		size_t getY() const;
		void setX(size_t);
		void setY(size_t);
};