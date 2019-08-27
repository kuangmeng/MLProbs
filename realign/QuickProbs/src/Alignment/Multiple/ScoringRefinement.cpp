#include "ScoringRefinement.h"
#include "EntropyEvaluator.h"

#include <algorithm>

using namespace std;
using namespace quickprobs;

#undef min
#undef max

bool quickprobs::ScoringRefinement::prepare(
	const float *seqsWeights, 
	const Array<float>& distances, 
	const Array<SparseMatrixType*> &sparseMatrices, 
	const ProbabilisticModel &model, 
	const MultiSequence &alignment)
{
	int numSeqs = alignment.count();
	EntropyEvaluator eval;

	columnScores.resize(alignment.length(), std::pair<int,float>(0,0));
	// iterate over columns
	for (int c = 0; c < columnScores.size(); ++c) {
		// iterate over sequences
		columnScores[c].first = c;
		columnScores[c].second = eval(alignment, c + 1);
	}

	std::stable_sort(columnScores.begin(), columnScores.end(), [](const std::pair<int,float> &a, const std::pair<int,float>& b){
		return a.second < b.second;
	});

	// erase elements
	int columnsUsed = (float)columnScores.size() * config->algorithm.refinement.columnFraction;

	int hi = std::min(std::max(columnsUsed, config->algorithm.refinement.iterations), (int)columnScores.size());

	if (hi > 0) {
		return true;
	}

	return false;
}
