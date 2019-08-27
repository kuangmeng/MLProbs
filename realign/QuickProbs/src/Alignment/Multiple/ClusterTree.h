#pragma once

#include "GuideTree.h"

#include "Common/Array.h"

namespace quickprobs 
{

class ClusterTree : public GuideTree
{
public:
	ClusterTree(Array<float>& distances);
	~ClusterTree();

protected:
	
	Array<float>& distances;
	
	//construct the cluster tree
	virtual void build();
};

};