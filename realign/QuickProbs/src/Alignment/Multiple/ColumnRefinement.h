#pragma once
#include "RefinementBase.h"
#include "Common/deterministic_random.h"

namespace quickprobs {

class ColumnRefinement : public RefinementBase
{
public:
	ColumnRefinement(
		std::shared_ptr<Configuration> config,
		std::shared_ptr<ConstructionStage> constructor) : RefinementBase(config, constructor) {
			if (config->algorithm.refinement.seed != 0) {
				engine.seed(config->algorithm.refinement.seed);
			}
	}

protected:	
	std::vector<float> posterior;
	
	std::mt19937 engine;

	std::vector<std::pair<int,float>> columnScores;

	virtual std::unique_ptr<MultiSequence> refine(
		const GuideTree& tree,
		const float *seqsWeights,
		const Array<float>& distances,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model, 
		std::unique_ptr<MultiSequence> alignment,
		int depth);

	virtual bool initialise(
		const float *seqsWeights,
		const Array<float>& distances,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model, 
		const MultiSequence &alignment);

	virtual bool finalise();

	virtual void split(
		const GuideTree& tree,
		const MultiSequence& alignment,
		std::set<int>& groupOne,
		std::set<int>& groupTwo);

	virtual void updateColumnScores(const MultiSequence &alignment, std::vector<std::pair<int,float>> &columnScores);
};

}