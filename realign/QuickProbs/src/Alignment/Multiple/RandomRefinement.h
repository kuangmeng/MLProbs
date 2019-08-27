#pragma once
#include "RefinementBase.h"

namespace quickprobs 
{

class RandomRefinement : public RefinementBase
{
public:
	RandomRefinement(
		std::shared_ptr<Configuration> config,
		std::shared_ptr<ConstructionStage> constructor) : RefinementBase(config, constructor) {}

protected:
	
	virtual void split(
		const GuideTree& tree,
		const MultiSequence& alignment,
		std::set<int>& groupOne,
		std::set<int>& groupTwo); 

	int genRandom(int m, int seed, bool init);
};

};