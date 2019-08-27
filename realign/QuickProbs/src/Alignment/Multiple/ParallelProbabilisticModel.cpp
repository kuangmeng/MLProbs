#include "Common/mathex.h"
#include "Common/Timer.h"
#include "DataStructures/Sequence.h"
#include "DataStructures/MultiSequence.h"
#include "ScoreType.h"
#include "ParallelProbabilisticModel.h"
#include "Configuration.h"

#include <omp.h>

#include <algorithm>
#include <memory>


#undef min

using namespace std;
using namespace quickprobs;

ParallelProbabilisticModel::ParallelProbabilisticModel(const Configuration &config) :
	ProbabilisticModel(config),
	numThreads(config.hardware.numThreads)
{
	localBufferSize = 512 * 512 * numThreads;
	localBuffer = new float[localBufferSize];
}

std::unique_ptr<std::vector<float>> quickprobs::ParallelProbabilisticModel::computeForwardMatrix(
	const Sequence& seq1, 
	const Sequence& seq2,
	float& total) const
{
	int height = seq1.GetLength() + 1;
	int width = seq2.GetLength() + 1;
	std::vector<float> *forwardPtr = new std::vector<float>(height * width, LOG_ZERO);
	this->computeForwardMatrix(seq1, seq2, forwardPtr->data(), total);
	return std::unique_ptr<std::vector<float>>(forwardPtr);
}

void quickprobs::ParallelProbabilisticModel::computeForwardMatrix(
	const Sequence& seq1, 
	const Sequence& seq2,
	float *forward,
	float& total) const
{
	int seq1Length = seq1.GetLength();
	int seq2Length = seq2.GetLength();
	int height = seq1Length + 1;
	int width = seq2Length + 1;
	auto data1 = seq1.getData();
	auto data2 = seq2.getData();

	// auxiliary matrices
	int layerOffset = width * 2;
	float* levels = new float[numStates * layerOffset];
	std::fill(levels, levels + numStates * layerOffset, LOG_ZERO);

	// initialization condition
	forward[1 * width + 1] = initialDistribution[0] + matchProb[(unsigned char) data1[1]][(unsigned char) data2[1]];

	for (int k = 0; k < numInsertStates; ++k) {
		levels[(2 * k + 1) * layerOffset + (1 * width + 0)] =
			initialDistribution[2 * k + 1] + insProb[(unsigned char) data1[1]][k];
		levels[(2 * k + 2) * layerOffset + (0 * width + 1)] =
			initialDistribution[2 * k + 2] + insProb[(unsigned char) data2[1]][k];
	}

	// remember offset for each index combination
	int ij = 0;
	int i1j = -width;
	int ij1 = -1;
	int i1j1 = -width - 1;

	int currentLevel = 0;
	int prevLevel = width;

	// compute forward scores
	for (int i = 0; i <= seq1Length; ++i) {
		unsigned char c1 = (unsigned char) data1[i];
		
		for (int j = 0; j <= seq2Length; ++j) {
			unsigned char c2 = (unsigned char) data2[j];

			int ij_level = currentLevel + j;
			int i1j_level = prevLevel + j;
			int ij1_level = currentLevel + j - 1;
			int i1j1_level = prevLevel + j - 1;

			if (i > 1 || j > 1) {
				if (i > 0 && j > 0) {
					forward[ij] = forward[i1j1] + transProb[0][0];
					for (int k = 1; k < numStates; ++k) {
						LOG_PLUS_EQUALS(forward[ij], levels[k * layerOffset + i1j1_level] + transProb[k][0]);
					}
					forward[ij] += matchProb[c1][c2];
				}
				if (i > 0) {
					for (int k = 0; k < numInsertStates; ++k) {
						int q = 2 * k + 1;
						levels[layerOffset * q + ij_level] = insProb[c1][k] + LOG_ADD(
							forward[i1j]+ transProb[0][q], levels[layerOffset * q + i1j_level] + transProb[q][q]);
					}
				}
				if (j > 0) {
					for (int k = 0; k < numInsertStates; ++k) {
						int q = 2 * k + 2;
						levels[layerOffset * q + ij_level] = insProb[c2][k] + LOG_ADD(
							forward[ij1] + transProb[0][q], levels[layerOffset * q + ij1_level] + transProb[q][q]);
					}
				}
			}

			++ij;
			++i1j;
			++ij1;
			++i1j1;
		}

		std::swap(currentLevel, prevLevel);
	}

	// get total forward probability
	total = LOG_ZERO; 
	LOG_PLUS_EQUALS(total, forward[width * height - 1] + initialDistribution[0]); // component from layer 0

	for (int k = 1; k < numStates; ++k) {
		LOG_PLUS_EQUALS(total, levels[layerOffset * k + prevLevel + seq2Length] + initialDistribution[k]); // components from other layers
	}

	delete [] levels;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> ParallelProbabilisticModel::computeBackwardMatrix(
	const Sequence & seq1, 
	const Sequence & seq2, 
	float& total) const
{
	const int height = seq1.GetLength() + 1;
	const int width = seq2.GetLength() + 1;
	std::vector<float> *backwardPtr = new std::vector<float>(height * width, LOG_ZERO);
	this->computeBackwardMatrix(seq1, seq2, backwardPtr->data(), total);
	return std::unique_ptr<std::vector<float>>(backwardPtr);
}


/// <summary>
/// See declaration for all the details.
/// </summary>
void ParallelProbabilisticModel::computeBackwardMatrix(
	const Sequence& seq1, 
	const Sequence& seq2,
	float *backward,
	float& total) const
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();
	int height = seq1Length + 1;
	int width = seq2Length + 1;
	auto data1 = seq1.getData();
	auto data2 = seq2.getData();

	// auxiliary matrices
	int layerOffset = width * 2;
	float* levels = new float[numStates * layerOffset];
	std::fill(levels, levels + numStates * layerOffset, LOG_ZERO);

	// initialization condition
	backward[width * height - 1] = initialDistribution[0];
	for (int k = 1; k < numStates; ++k)
		levels[layerOffset * k + seq2Length] = initialDistribution[k];

	// remember offset for each index combination
	int ij = height * width - 1;
	int i1j1 = ij + width + 1;

	int currentLevel = 0;
	int nextLevel = width;

	// compute backward scores
	for (int i = seq1Length; i >= 0; --i) {
		unsigned char c1 = (i == seq1Length) ? '~' : (unsigned char) data1[i + 1];
		for (int j = seq2Length; j >= 0; --j) {
			unsigned char c2 = (j == seq2Length) ? '~' : (unsigned char) data2[j + 1];

			int ij_level = currentLevel + j;
			int i1j_level = nextLevel + j;
			int ij1_level = currentLevel + j + 1;
		
			// reset current cell
			if (i < seq1Length || j < seq2Length) {
				for (int k = 1; k < numStates; ++k)
					levels[layerOffset * k + ij_level] = LOG_ZERO; 
			}

			if (i < seq1Length && j < seq2Length) {
				const float ProbXY = backward[i1j1] + matchProb[c1][c2];
				LOG_PLUS_EQUALS(backward[ij], ProbXY + transProb[0][0]);
				for (int k = 1; k < numStates; ++k)
					LOG_PLUS_EQUALS(levels[layerOffset * k + ij_level], ProbXY + transProb[k][0]);
			}
			if (i < seq1Length) {
				for (int k = 0; k < numInsertStates; ++k) {
					int q = 2 * k + 1;
					LOG_PLUS_EQUALS(backward[ij], levels[layerOffset * q + i1j_level] + insProb[c1][k] + transProb[0][q]);
					LOG_PLUS_EQUALS(levels[layerOffset * q + ij_level], levels[layerOffset * q + i1j_level] + insProb[c1][k] + transProb[q][q]);
				}
			}
			if (j < seq2Length) {
				for (int k = 0; k < numInsertStates; ++k) {
					int q = 2 * k + 2;
					LOG_PLUS_EQUALS(backward[ij], levels[layerOffset * q + ij1_level] + insProb[c2][k] + transProb[0][q]);
					LOG_PLUS_EQUALS(levels[layerOffset * q + ij_level], levels[layerOffset * q + ij1_level] + insProb[c2][k] + transProb[q][q]);
				}
			}

			--ij;
			--i1j1;
		}

		std::swap(currentLevel, nextLevel);
	}

	total = initialDistribution[0] + matchProb[data1[1]][data2[1]] + backward[width + 1]; // get element from [1,1]

	for (int k = 0; k < numInsertStates; ++k) {
		LOG_PLUS_EQUALS(total, initialDistribution[2*k+1] + insProb[data1[1]][k] + levels[layerOffset * (2*k+1) + currentLevel]);
		LOG_PLUS_EQUALS(total, initialDistribution[2*k+2] + insProb[data2[1]][k] + levels[layerOffset * (2*k+2) + nextLevel + 1]);
	}

	delete [] levels;
}


/// <summary>
/// See declaration for all the details.
/// </summary>
void ParallelProbabilisticModel::computePosteriorMatrix(
	const Sequence & seq1, 
	const Sequence & seq2,
	const float* forward, 
	const float* backward,
	float totalProb,
	float *posterior) const
{
	const int seq1Length = seq1.GetLength();
	const int seq2Length = seq2.GetLength();

	int ij = 0;
	if (totalProb == 0) {
		totalProb = 1.0f;
	}
	
	float* out = posterior;
	const float* m1 = forward;
	const float* m2 = backward;

	int layerSize = (seq1Length + 1) * (seq2Length + 1);
	for (int id = 0; id < layerSize; ++id) {
		*out = EXP(std::min(LOG_ONE, *m1 + *m2 - totalProb));
		++out;
		++m1;
		++m2;
	}

	posterior[0] = 0;
}


/// <summary>
/// See declaration for all the details.
/// </summary>
std::unique_ptr<std::vector<float>> ParallelProbabilisticModel::buildPosterior(
	const float* seqsWeights, 
	const Array<float>& distances,
	const quickprobs::MultiSequence & align1, 
	const quickprobs::MultiSequence & align2,
	const Array<SparseMatrixType*> &sparseMatrices,
	float selectivity) const 
{
	// OpenMP version
	const int seq1Length = align1.GetSequence(0)->GetLength();
	const int seq2Length = align2.GetSequence(0)->GetLength();

	::size_t numElements = (seq1Length + 1) * (seq2Length + 1);

	// rounded size to enable loop unrolling
	std::vector<float>* posterior = new std::vector<float>(numElements, 0);

	this->buildPosterior(
		seqsWeights, distances, align1, align2, sparseMatrices, selectivity, posterior->data());

	return std::unique_ptr<std::vector<float>>(posterior);
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void ParallelProbabilisticModel::buildPosterior(
	const float* seqsWeights, 
	const Array<float>& distances,
	const quickprobs::MultiSequence & align1, 
	const quickprobs::MultiSequence & align2,
	const Array<SparseMatrixType*> &sparseMatrices,
	float selectivity,
	float * posterior) const 
{
	TIMER_CREATE(timer);
	
	TIMER_START(timer);
	int numPairs = align1.count() * align2.count();
	const int seq1Length = align1.length();
    const int seq2Length = align2.length();
	int * seqsMappings1 = new int[align1.count() * (seq1Length + 1)];
	int * seqsMappings2 = new int[align2.count() * (seq2Length + 1)];

	// rounded size to enable loop unrolling
	::size_t numElements = (seq1Length + 1) * (seq2Length + 1);
	std::memset(posterior, 0, numElements * sizeof(float));

	//compute the total sum of all weights
	double totalWeights = 0;
    for (int i = 0; i < align1.count(); ++i){
		int first = align1.GetSequence(i)->GetLabel();
	  	double w1 = seqsWeights[first];
      	for (int j = 0; j < align2.count(); ++j){
 			int second = align2.GetSequence(j)->GetLabel();
			
			if (distances[first][second] <= selectivity) {
				double w2 = seqsWeights[second];
				totalWeights += w1 * w2;
			}
	  	}
	}

	//#pragma omp parallel for default(shared) schedule(dynamic)
	for(int i = 0; i < align1.count(); ++i) {
		align1.GetSequence(i)->getMapping(seqsMappings1 + i * (seq1Length + 1));
	}
	//#pragma omp parallel for default(shared) schedule(dynamic)
	for(int j = 0; j < align2.count(); ++j) {
		align2.GetSequence(j)->getMapping(seqsMappings2 + j * (seq2Length + 1));
	}

	TIMER_STOP(timer);
	timePreparation += timer.seconds();


	TIMER_START(timer);
	// Serial variant
	if (numThreads == 1 || numPairs < numThreads * 2) {
		for (int i = 0; i < align1.count(); i++) {
			int first = align1.GetSequence(i)->GetLabel();
			const int* mapping1 = seqsMappings1 + i * (seq1Length + 1);
			double w1 = seqsWeights[first];

			// for each t in align2
			for (int j = 0; j < align2.count(); j++) {
				int second = align2.GetSequence(j)->GetLabel();
				const int* mapping2 = seqsMappings2 + j * (seq2Length + 1);
				double w2 = seqsWeights[second];	

				float w = (float)((w1 * w2) / totalWeights);
				const auto matrix = sparseMatrices[first][second];
				int currentSeq1Length = matrix->getSeq1Length();

				for (int ii = 1; ii <= currentSeq1Length; ++ii) {
					auto row = matrix->getRowPtr(ii);
					auto rowEnd = row + matrix->getRowSize(ii);
					ptrdiff_t denseBase = (mapping1)[ii] * (seq2Length + 1);

					// add in all relevant values
					for (; row < rowEnd; ++row) {
						float v = row->getValue();
						int col = row->getColumn();
						int id = denseBase + (mapping2)[col];
						posterior[id] += w * v;
					}
				}
			}
		}
	} else {
		// parallel variant
		#pragma omp parallel default(shared)
		{
			int tid = omp_get_thread_num();
			
			for (int i = 0; i < align1.count(); i++) {
				int first = align1.GetSequence(i)->GetLabel();
				const int* mapping1 = seqsMappings1 + i * (seq1Length + 1);
				double w1 = seqsWeights[first];

				// fixme: just to take length of ungapped sequence
				int second = align2.GetSequence(0)->GetLabel();
				int currentSeq1Length = sparseMatrices[first][second]->getSeq1Length();

				int blockSize = mathex::ceildiv(currentSeq1Length, numThreads);
				int blockBegin = 1 + blockSize * tid;
				int blockEnd = min(currentSeq1Length + 1, blockBegin + blockSize);

				// for each t in align2
				for (int j = 0; j < align2.count(); j++) {
					int second = align2.GetSequence(j)->GetLabel();
					const int* mapping2 = seqsMappings2 + j * (seq2Length + 1);
					double w2 = seqsWeights[second];	
				
					float w = (float)((w1 * w2) / totalWeights);
					const auto matrix = sparseMatrices[first][second];
					
					
					for (int ii = blockBegin; ii < blockEnd; ++ii) {

					//for (int ii = 1; ii <= currentSeq1Length; ++ii) {
					//	if (mapping1[ii] % numThreads == tid) {
							auto row = matrix->getRowPtr(ii);
							auto rowEnd = row + matrix->getRowSize(ii);
							ptrdiff_t denseBase = (mapping1)[ii] * (seq2Length + 1);
						
							// add in all relevant values
							for (; row < rowEnd; ++row) {
								float v = row->getValue();
								int col = row->getColumn();
								int id = denseBase + (mapping2)[col];
								posterior[id] += w * v;
							}
						//}
					}
				}

				#pragma omp barrier
			}
		}
	}

	TIMER_STOP(timer);
	timeCalculaton += timer.seconds();

	TIMER_START(timer);
	delete [] seqsMappings1;
	delete [] seqsMappings2;
	TIMER_STOP(timer);
	timePreparation += timer.seconds();
}