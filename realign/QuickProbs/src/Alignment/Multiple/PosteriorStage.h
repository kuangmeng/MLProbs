#pragma once
#include "IAlgorithmStage.h"
#include "BufferSet.h"
#include "Alignment/DataStructures/SparseMatrixType.h"

#include "Common/Array.h"

namespace quickprobs
{
// forward declare
class ISequenceSet;
class ProbabilisticModel;
class PartitionFunction;
class Configuration;
class Sequence;
class MultiSequence;

class PosteriorStage : public IAlgorithmStage
{
public:
	const std::shared_ptr<ProbabilisticModel> getModel() const { return model; } 

	std::shared_ptr<ProbabilisticModel> getModel() { return model; } 

	PosteriorStage(std::shared_ptr<Configuration> config);

	virtual void operator()(
		ISequenceSet& set, 
		Array<float>& distances,
		Array<SparseMatrixType*>& matrices);

	virtual void computePairwise(
		const Sequence& seq1, 
		const Sequence& seq2,
		BufferSet& buffers, 
		float& distance);

protected:
	std::shared_ptr<ProbabilisticModel> model;

	std::shared_ptr<PartitionFunction> function;

	size_t sparseBytes;

	size_t denseBytes;
	
	virtual void run(
		ISequenceSet& set, 
		Array<float>& distances,
		Array<SparseMatrixType*>& matrices);

	float combineMatrices(
		int seq1Length,
		int seq2Length,
		const float* input1, 
		const float* input2, 
		float* out);
};

};