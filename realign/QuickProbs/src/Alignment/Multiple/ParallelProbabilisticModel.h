#pragma once
#include <memory>
#include "ProbabilisticModel.h"
#include "ProbabilisticParams.h"
#include "AuxiliaryStructures.h"

namespace quickprobs
{

/// <summary>
/// This class calculates all the probabilistic stuff using kernels. It is used basically for
/// testing kernel implementations.
/// </summary>
class ParallelProbabilisticModel : public quickprobs::ProbabilisticModel
{
public:
	virtual void setNumThreads(int numThreads) { this->numThreads = numThreads; }

	virtual int getNumThreads() const { return numThreads; }

	ParallelProbabilisticModel(const Configuration &config);

	virtual ~ParallelProbabilisticModel() { delete [] localBuffer; }

	virtual	std::unique_ptr<std::vector<float>> computeForwardMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		float& total) const; 

	virtual void computeForwardMatrix(
		const Sequence& seq1, 
		const Sequence& seq2,
		float *forward,
		float& total) const;

	virtual	std::unique_ptr<std::vector<float>> computeBackwardMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		float& total) const;

	virtual void computeBackwardMatrix(
		const Sequence& seq1, 
		const Sequence& seq2,
		float *backward,
		float& total) const;

	virtual void computePosteriorMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		const float* forward, 
		const float* backward,
		float totalProb,
		float *posterior) const;

	virtual std::unique_ptr<std::vector<float>> buildPosterior(
		const float* seqsWeights,
		const Array<float>& distances,
		const MultiSequence& align1, 
		const MultiSequence& align2,
		const Array<SparseMatrixType*> &sparseMatrices,
		float selectivity) const;

	virtual void buildPosterior(
		const float* seqsWeights,
		const Array<float>& distances,
		const MultiSequence& align1, 
		const MultiSequence& align2,
		const Array<SparseMatrixType*> &sparseMatrices,
		float selectivity,
		float * posterior) const;

protected:
	int numThreads;

	mutable float* localBuffer;

	mutable ::size_t localBufferSize;
	
};

};