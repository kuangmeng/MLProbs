#pragma once
#include <limits>
#undef max
#undef min

namespace quickprobs 
{

class ITasksWave
{
public:
	std::map<int, int> verticalLengths;
	std::map<int, int> horizontalLengths;

	virtual ::size_t size() const = 0;
	
	int maxLength() const { return std::max(maxVertical(), maxHorizontal()); }

	int minVertical() const { return verticalLengths.cbegin()->first; }
	int maxVertical() const { return verticalLengths.crbegin()->first; }

	int minHorizontal() const { return horizontalLengths.cbegin()->first; }
	int maxHorizontal() const { return horizontalLengths.crbegin()->first; }
	

protected:
	
};
};