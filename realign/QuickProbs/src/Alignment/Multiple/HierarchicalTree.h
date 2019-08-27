#pragma once
#include "GuideTree.h"

namespace quickprobs 
{


class HierarchicalTree : public GuideTree
{
public:
	HierarchicalTree(std::vector<std::vector<float>>& distances);
	~HierarchicalTree();

protected:
	std::vector<std::vector<float>>& distances;

	virtual void build();
};

}