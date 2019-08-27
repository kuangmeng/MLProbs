#pragma once
#include "GuideTree.h"

namespace quickprobs {

class SLinkTree : public GuideTree
{
public:
	SLinkTree(Array<float>& distances);

protected:
	Array<float>& distances;

	virtual void build();
};

}