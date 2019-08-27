#include "Common/mathex.h"
#include "Common/Log.h"

#include "DataStructures/Sequence.h"
#include "DataStructures/MultiSequence.h"
#include "DataStructures/SparseHelper.h"
#include "Pairwise/ProteinHmm5.h"
#include "Pairwise/AminoAcidMatrices.h"

#include "ParallelProbabilisticModel.h"
#include "PartitionFunction.h"
#include "PosteriorStage.h"
#include "ProbabilisticParams.h"
#include "ScoreType.h"


#include <omp.h>

#undef min
#undef max

using namespace std;
using namespace quickprobs;
using namespace mathex;

/// <summary>
/// See declaration for all the details.
/// </summary>
PosteriorStage::PosteriorStage(std::shared_ptr<Configuration> config) : IAlgorithmStage(config)
{
	// partition function configuration
	function = std::shared_ptr<PartitionFunction>(new PartitionFunction(*config));

	// configure hmms
	model = std::shared_ptr<quickprobs::ProbabilisticModel>(new quickprobs::ParallelProbabilisticModel(*config));
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void PosteriorStage::operator()(
	ISequenceSet& set,
	Array<float>& distances,
	Array<SparseMatrixType*>& sparseMatrices)
{
	TIMER_CREATE(timer);
	TIMER_START(timer);
	run(set, distances, sparseMatrices);
	TIMER_STOP(timer);

	size_t elemsCount = config->io.enableVerbose ? SparseHelper::totalElements(sparseMatrices) : 0;
	double elemsSum = config->io.enableVerbose ? SparseHelper::sumOfElements(sparseMatrices) : 0;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void PosteriorStage::run(
	ISequenceSet& set,
	Array<float>& distances,
	Array<SparseMatrixType*>& sparseMatrices)
{
	const int numSeqs = set.count();

	//calculate sequence pairs for openmp model
	omp_set_num_threads(config->hardware.numThreads);
	int pairIdx = 0;
	int numPairs = (numSeqs - 1) * numSeqs / 2;
	std::vector<std::pair<int,int>> seqsPairs(numPairs);
	::size_t maxSeqLength = 0;

	for(int a = 0; a < numSeqs; a++){
		if (set.GetSequence(a)->GetLength() > maxSeqLength) {
			maxSeqLength = set.GetSequence(a)->GetLength();
		}

		for(int b = a + 1; b < numSeqs; b++){
			seqsPairs[pairIdx].first = a;
			seqsPairs[pairIdx].second = b;
			pairIdx++;
		}
	}

	int maxLayerSize = (maxSeqLength + 1) * (maxSeqLength + 1);
	std::vector<std::shared_ptr<BufferSet>> bufferSets(config->hardware.numThreads);
	for (auto& set : bufferSets) {
		set = std::shared_ptr<BufferSet>(new BufferSet(maxLayerSize));
	}


	// do all pairwise alignments for posterior probability matrices
	size_t denseBytes = 0;
	size_t sparseBytes = 0;
	#pragma omp parallel for default(shared) schedule(dynamic) reduction(+: denseBytes, sparseBytes)
	for(int pairIdx = 0; pairIdx < numPairs; pairIdx++){
		int i = seqsPairs[pairIdx].first;
		int j = seqsPairs[pairIdx].second;

		int tid = omp_get_thread_num();
		BufferSet& bufset = *bufferSets[tid];

		Sequence* seq1 = set.GetSequence(i);
		Sequence* seq2 = set.GetSequence(j);

		computePairwise(*seq1, *seq2, bufset, distances[i][j]);
		denseBytes += maxLayerSize * sizeof(float);

		//// compute sparse representations
		auto m = new SparseMatrixType(seq1->GetLength(), seq2->GetLength(), bufset.f0(), config->algorithm.posteriorCutoff);
		sparseMatrices[i][j] = m;
		sparseMatrices[j][i] = m->computeTranspose();
		distances[j][i] = distances[i][j];
		sparseBytes += sparseMatrices[i][j]->getNumCells() * sizeof(SparseMatrixType::cell_type) +
			sparseMatrices[j][i]->getHeight() * sizeof(int);
	}

}

/// <summary>
/// See declaration for all the details.
/// </summary>

void quickprobs::PosteriorStage::computePairwise(const Sequence& seq1, const Sequence& seq2, BufferSet& buffers, float& distance)
{
	//compute posterior probability matrix from partition function
	int layerSize = (seq1.GetLength() + 1) * (seq2.GetLength() + 1);
	std::fill_n(buffers.f2(), layerSize, 0);
	std::fill_n(buffers.d01(), layerSize, 0);

	distance = 0;
	function->computePosteriorMatrix(seq1, seq2, buffers.d01(), buffers.f2());

	// compute forward and backward probabilities
	float totalForward = 0;
	float totalBackward = 0;

	std::fill_n(buffers.f0(), layerSize, LOG_ZERO);
	std::fill_n(buffers.f1(), layerSize, LOG_ZERO);

	model->computeForwardMatrix(seq1, seq2, buffers.f0(), totalForward);
	model->computeBackwardMatrix(seq1, seq2, buffers.f1(), totalBackward);
	float total = (totalForward + totalBackward) / 2;

	// compute posterior probability matrix from HMM
	model->computePosteriorMatrix(seq1, seq2, buffers.f0(), buffers.f1(), total, buffers.f1());

	//merge the two posterior matrices
	distance = combineMatrices(
		seq1.GetLength(),
		seq2.GetLength(),
		buffers.f1(),
		buffers.f2(),
		buffers.f0());
}

float quickprobs::PosteriorStage::combineMatrices(
	int seq1Length,
	int seq2Length,
	const float* input1,
	const float* input2,
	float* out)
{
	float *twoRows = new float[(seq2Length + 1) * 2];
	float *oldRow = twoRows;
	float *newRow = twoRows + seq2Length + 1;

	for(int i = 0; i <= seq1Length; ++i) {
		for(int j = 0; j <= seq2Length; ++j) {
			// fixed:
			if (i == 0 || j == 0) {
				*out = 0;
				newRow[j] = 0;
			} else {
				float v1 = *input1;
				float v2 = *input2;
				*out = sqrt((v1 * v1 + v2 * v2) * 0.5f);
				newRow[j] = mathex::max(*out + oldRow[j-1], newRow[j-1], oldRow[j]);
			}

			++input1;
			++input2;
			++out;
		}

		// swap rows
		float *temp = oldRow;
		oldRow = newRow;
		newRow = temp;
	}

	float total = oldRow[seq2Length];
	delete [] twoRows;

	float distance = 1.0f - total / std::min(seq1Length, seq2Length);
	return distance;
}
