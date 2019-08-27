#pragma once
#include "GuideTree.h"

namespace quickprobs 
{
class MultiSequence;

class NewickTree : public GuideTree
{
public:
	NewickTree(MultiSequence& sequences, std::string description);

protected:
	
	MultiSequence &sequences;
	
	std::string description;

	virtual void build();
};

};