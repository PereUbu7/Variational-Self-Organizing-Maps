class SomIndex
{
	protected:
		int x,y;
	
	public:
		SomIndex(int, int);
		~SomIndex() = default;
		//int getSomIndex(Som);
		unsigned int getX() const;
		unsigned int getY() const;
		void setX(int);
		void setY(int);
};