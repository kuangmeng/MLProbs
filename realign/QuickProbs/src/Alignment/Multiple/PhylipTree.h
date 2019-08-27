#pragma once
#include "NewickTree.h"

namespace quickprobs 
{
class MultiSequence;

enum PhylogenyType {
	LIKELIHOOD = 0,
	PARSIMONY
};

class PhylipTree : public NewickTree
{
public:
	PhylogenyType type;
	
	PhylipTree(MultiSequence& alignment, PhylogenyType type);

protected:
	virtual void build();

};

};
