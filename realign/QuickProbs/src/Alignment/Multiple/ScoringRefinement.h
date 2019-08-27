#pragma once
#include "ColumnRefinement.h"

namespace quickprobs {

class ScoringRefinement : public ColumnRefinement
{
public:
	ScoringRefinement(
		std::shared_ptr<Configuration> config,
		std::shared_ptr<ConstructionStage> constructor) : ColumnRefinement(config, constructor) {}

protected:
	virtual bool prepare(
		const float *seqsWeights,
		const Array<float>& distances,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model, 
		const MultiSequence &alignment);
};

}