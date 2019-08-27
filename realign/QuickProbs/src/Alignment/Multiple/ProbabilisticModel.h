/////////////////////////////////////////////////////////////////
// ProbabilisticModel.h
//
// Routines for (1) posterior probability computations
//              (2) chained anchoring
//              (3) maximum weight trace alignment
/////////////////////////////////////////////////////////////////

// Adam Gudys modifications:
// 06.03.2012: Some methods made virtual.
// 11.04.2012: SymbolCount constant added.
#pragma once
#include <memory>
#include <vector>

#include "Alignment/DataStructures/SparseMatrixType.h"
#include "Alignment/Multiple/Common.h"
#include "ProbabilisticParams.h"

#include "Common/Array.h"

namespace quickprobs
{
// forward declarations
class Sequence;
class MultiSequence;
class Configuration;

/////////////////////////////////////////////////////////////////
// ProbabilisticModel
//
// Class for storing the parameters of a probabilistic model and
// performing different computations based on those parameters.
// In particular, this class handles the computation of
// posterior probabilities that may be used in alignment.
/////////////////////////////////////////////////////////////////

class ProbabilisticModel {
protected:	
	int numStates;
	int numInsertStates;
	
	float initialDistribution[5];	// holds the initial probabilities for each state
	float transProb[5][5];			// holds all state-to-state transition probabilities
	float matchProb[256][256];     // emission probabilities for match states
	float insProb[256][5];			// emission probabilities for insert states
	
public:
	ProbabilisticParams params;
	
	mutable double timeExecution;
	mutable double timePreparation;
	mutable double timeCalculaton;
	mutable double timeReduction;

	virtual void setNumThreads(int numThreads) { throw std::runtime_error("Not implemented!"); }

	virtual int getNumThreads() const { return 1; }

	float getMatchProb(int i, int j) const { return matchProb[i][j]; }

	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ProbabilisticModel()
	//
	// Constructor.  Builds a new probabilistic model using the
	// given parameters.
	/////////////////////////////////////////////////////////////////
	ProbabilisticModel(const Configuration& config);

	virtual ~ProbabilisticModel() {}
	
	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputeForwardMatrix()
	//
	// Computes a set of forward probability matrices for aligning
	// seq1 and seq2.
	//
	// For efficiency reasons, a single-dimensional floating-point
	// array is used here, with the following indexing scheme:
	//
	//    forward[i + NumMatrixTypes * (j * (seq2Length+1) + k)]
	//    refers to the probability of aligning through j characters
	//    of the first sequence, k characters of the second sequence,
	//    and ending in state i.
	/////////////////////////////////////////////////////////////////
	virtual	std::unique_ptr<std::vector<float>> computeForwardMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		float& total) const; 

	virtual void computeForwardMatrix(
		const Sequence& seq1, 
		const Sequence& seq2,
		float *forward,
		float& total) const;
	
	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputeBackwardMatrix()
	//
	// Computes a set of backward probability matrices for aligning
	// seq1 and seq2.
	//
	// For efficiency reasons, a single-dimensional floating-point
	// array is used here, with the following indexing scheme:
	//
	//    backward[i + NumMatrixTypes * (j * (seq2Length+1) + k)]
	//    refers to the probability of starting in state i and
	//    aligning from character j+1 to the end of the first
	//    sequence and from character k+1 to the end of the second
	//    sequence.
	/////////////////////////////////////////////////////////////////
	virtual	std::unique_ptr<std::vector<float>> computeBackwardMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		float& total) const;

	virtual void computeBackwardMatrix(
		const Sequence& seq1, 
		const Sequence& seq2,
		float *backward,
		float& total) const;

	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputeTotalProbability()
	//
	// Computes the total probability of an alignment given
	// the forward and backward matrices.
	/////////////////////////////////////////////////////////////////
	virtual	float computeTotalProbability(
		int seq1Length, 
		int seq2Length,
		const float* forward, 
		const float* backward) const;


	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputePosteriorMatrix()
	//
	// Computes the posterior probability matrix based on
	// the forward and backward matrices.
	/////////////////////////////////////////////////////////////////
	virtual	std::unique_ptr<std::vector<float>> computePosteriorMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		const std::vector<float> &forward, 
		const std::vector<float> &backward,
		float totalProb) const;

	virtual void computePosteriorMatrix(
		const Sequence & seq1, 
		const Sequence & seq2,
		const float* forward, 
		const float* backward,
		float totalProb,
		float* posterior) const;

	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputeAlignment()
	//
	// Computes an alignment based on the given posterior matrix.
	// This is done by finding the maximum summing path (or
	// maximum weight trace) through the posterior matrix.  The
	// final alignment is returned as a pair consisting of:
	//    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
	//        denote insertions in one of the two sequences and
	//        B's denote that both sequences are present (i.e.
	//        matches).
	//    (2) a float indicating the sum achieved
	/////////////////////////////////////////////////////////////////
	virtual std::pair<std::vector<char>*, float> computeAlignment(
		int seq1Length,
		int seq2Length, 
	    float* posterior) const;

	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputeAlignmentWithGapPenalties()
	//
	// Similar to ComputeAlignment() except with gap penalties.
	/////////////////////////////////////////////////////////////////
	virtual std::pair<std::vector<char> *, float> computeAlignmentWithGapPenalties(
			const MultiSequence & align1, 
			const MultiSequence & align2, 
			const std::vector<float> &posterior,
			int numSeqs1, 
			int numSeqs2, 
			float gapOpenPenalty,
			float gapContinuePenalty) const;

	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::ComputeViterbiAlignment()
	//
	// Computes the highest probability pairwise alignment using the
	// probabilistic model.  The final alignment is returned as a
	//  pair consisting of:
	//    (1) a string (e.g., XXXBBXXXBBBBBBYYYYBBB) where X's and
	//        denote insertions in one of the two sequences and
	//        B's denote that both sequences are present (i.e.
	//        matches).
	//    (2) a float containing the log probability of the best
	//        alignment (not used)
	/////////////////////////////////////////////////////////////////
	virtual std::pair<std::vector<char> *, float> computeViterbiAlignment(
		const Sequence & seq1, 
		const Sequence & seq2) const;

	/////////////////////////////////////////////////////////////////
	// ProbabilisticModel::BuildPosterior()
	//
	// Builds a posterior probability matrix needed to align a pair
	// of alignments.  Mathematically, the returned matrix M is
	// defined as follows:
	//    M[i,j] =     sum          sum      f(s,t,i,j)
	//             s in align1  t in align2
	// where
	//                  [  P(s[i'] <--> t[j'])
	//                  [       if s[i'] is a letter in the ith column of align1 and
	//                  [          t[j'] it a letter in the jth column of align2
	//    f(s,t,i,j) =  [
	//                  [  0    otherwise
	//
	/////////////////////////////////////////////////////////////////
	virtual std::unique_ptr<std::vector<float>> buildPosterior(
		const MultiSequence & align1, 
		const MultiSequence & align2,
		const Array<SparseMatrixType*> &sparseMatrices,
		float cutoff) const;

	virtual std::unique_ptr<std::vector<float>> buildPosterior(
		const float* seqsWeights, 
		const MultiSequence & align1,
		const MultiSequence & align2,
		const Array<SparseMatrixType*> &sparseMatrices,
		float cutoff) const; 
	};
};
