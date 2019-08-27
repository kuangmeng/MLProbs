#pragma once
#include "IAlgorithmStage.h"
#include "Selectivity.h"
#include "Alignment/DataStructures/SparseMatrixType.h"

#include "Common/Array.h"

namespace quickprobs 
{

// Forward declarations
class MultiSequence;
class ISequenceSet;

// Class declaration
class ConsistencyStage : public IAlgorithmStage
{
public:
	ConsistencyStage(std::shared_ptr<Configuration> config);
	
	virtual void operator()(
		const float* seqsWeights,
		ISequenceSet& set,
		Array<float>& distances,
		Array<SparseMatrixType*>& sparseMatrices);

protected:
	Selectivity selectivity;

	int iterations;
	float selfweight;
	
	virtual void run(
		const float* seqsWeights,
		ISequenceSet& set, 
		Array<float>& distances,
		Array<SparseMatrixType*>& sparseMatrices);

	virtual Array<SparseMatrixType*> doRelaxation(
		const float* seqsWeights, 
		const ISequenceSet *sequences, 
		const Array<float>& distances,
		const Array<SparseMatrixType*> &sparseMatrices,
		bool filterFlag);

	virtual void relax(float weight, const SparseMatrixType *matXZ, const SparseMatrixType *matZY, std::vector<float>& posterior);
	virtual void relaxTransposed(float weight, const SparseMatrixType *matZX, const SparseMatrixType *matZY, std::vector<float> &posterior);
};

};