#pragma once 
#include "RefinementBase.h"

namespace quickprobs {

class TreeRefinement : public RefinementBase {
public:
	TreeRefinement(
			std::shared_ptr<Configuration> config,
			std::shared_ptr<ConstructionStage> constructor) : RefinementBase(config, constructor) {}

protected:
	std::mt19937 engine;
	
	virtual void split(
		const GuideTree& tree,
		const MultiSequence& alignment,
		std::set<int>& groupOne,
		std::set<int>& groupTwo);

};

};