#pragma once
#include "IAlgorithmStage.h"
#include "GuideTree.h"
#include "ProbabilisticModel.h"
#include "Alignment/DataStructures/MultiSequence.h"
#include "Alignment/DataStructures/SparseMatrixType.h"

namespace quickprobs 
{

class ConstructionStage : public IAlgorithmStage
{
public:
	// two timers for experimental purposes
	double timeMatrixAddition;
	double timeDynamicProgramming;
	
	ConstructionStage(std::shared_ptr<Configuration> config) 
		: IAlgorithmStage(config), selectivity(config->algorithm.finalSelectivity) {}
	
	virtual std::unique_ptr<quickprobs::MultiSequence> operator()(
		const float *seqsWeights,
		const Array<float>& distances,
		GuideTree* tree, 
		const MultiSequence& sequences,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model);

	virtual std::unique_ptr<quickprobs::MultiSequence> alignAlignments(
		const float *seqsWeights,
		const Array<float>& distances,
		const MultiSequence& align1, 
		const MultiSequence& align2,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model,
		float * buffer);

protected:
	
	float selectivity;

	virtual std::unique_ptr<quickprobs::MultiSequence> run(
		const float *seqsWeights,
		const Array<float>& distances,
		GuideTree* tree, 
		const MultiSequence& sequences,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model);

	virtual std::unique_ptr<quickprobs::MultiSequence> processTree(
		const float *seqsWeights,
		const Array<float>& distances,
		const Node& tree, 
		const MultiSequence& sequences,
		const Array<SparseMatrixType*> &sparseMatrices,
		const ProbabilisticModel &model);
};

};