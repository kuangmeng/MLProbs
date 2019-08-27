#include "ColumnRefinement.h"
#include "DataStructures/Sequence.h"
#include "DataStructures/MultiSequence.h"

#include "Common/mathex.h"
#include "Common/Log.h"

#include <random>
#include <utility>
#include <cmath>

#undef min
#undef max

using namespace quickprobs;


std::unique_ptr<MultiSequence> quickprobs::ColumnRefinement::refine(
	const GuideTree& tree,
	const float *seqsWeights, 
	const Array<float>& distances, 
	const Array<SparseMatrixType*> &sparseMatrices, 
	const ProbabilisticModel &model, 
	std::unique_ptr<MultiSequence> alignment, 
	int depth)
{
	std::set<int> groupOne, groupTwo;
	split(tree, *alignment, groupOne, groupTwo);

	if (groupOne.size() > 0 && groupTwo.size() > 0) {
		auto profileOne = alignment->extractSubset(groupOne); 
		auto profileTwo = alignment->extractSubset(groupTwo);
		
		// refine recursively
		int maxDepth = std::min(config->algorithm.refinement.maxDepth, (int)std::log2(distances.size()));
		if (depth < maxDepth) {
			profileOne = refine(tree, seqsWeights, distances, sparseMatrices, model, std::move(profileOne), depth + 1);
			profileTwo = refine(tree, seqsWeights, distances, sparseMatrices, model, std::move(profileTwo), depth + 1);
		}

		// resize posterior buffer when necessary
		::size_t memNeeded = (profileOne->length() + 1) * (profileTwo->length() + 1);
		if (posterior.size() < memNeeded) {
			posterior.resize(std::max(posterior.size() * 2, memNeeded));
		}

		auto candidate = constructor->alignAlignments(seqsWeights, distances, *profileOne, *profileTwo, sparseMatrices, model, posterior.data());

		if (checkAcceptance(*alignment, *candidate)) {
			alignment = std::move(candidate);
		} 
	} 

	return std::move(alignment);
}



bool quickprobs::ColumnRefinement::initialise(
	const float *seqsWeights,
	const Array<float>& distances,
	const Array<SparseMatrixType*> &sparseMatrices,
	const ProbabilisticModel &model, 
	const MultiSequence &alignment)
{
	int numSeqs = alignment.count();
	posterior.resize((alignment.length() + 1) * (alignment.length() + 1));

	updateColumnScores(alignment, columnScores);

	int columnsUsed = (float)columnScores.size() * std::abs(config->algorithm.refinement.columnFraction);


	int hi = std::min(std::max(columnsUsed, config->algorithm.refinement.iterations), (int)columnScores.size());
	return hi > 0;
}


bool quickprobs::ColumnRefinement::finalise() {
	std::size_t h = 0;
	double sum = 0.0;
	for (float e : posterior) {
		h = h ^ std::hash<float>()(e);
		sum += e;
	}
	
	return true;
}

void quickprobs::ColumnRefinement::split(
	const GuideTree& tree,
	const MultiSequence& alignment,
	std::set<int>& groupOne,
	std::set<int>& groupTwo)
{
	int numSeqs = alignment.count();
	updateColumnScores(alignment, columnScores);

	int columnsUsed = (float)columnScores.size() * std::abs(config->algorithm.refinement.columnFraction);
	int lo, hi;

	if (config->algorithm.refinement.columnFraction > 0) {
		lo = 0;
		hi = std::min(std::max(columnsUsed, config->algorithm.refinement.iterations), (int)columnScores.size());
	} else {
		lo = std::max((size_t)0, columnScores.size() - std::max(columnsUsed, config->algorithm.refinement.iterations));
		hi = columnScores.size();
	}
	
	if (hi > 0) {
		det_uniform_int_distribution<int> distribution(lo, hi - 1);
	
		int rnd = distribution(this->engine);
		int divisionColumn = std::min((size_t)this->columnScores[rnd].first, alignment.length() - 1);

		for (int i = 0; i < numSeqs; ++i) {
			if (alignment.GetSequence(i)->GetPosition(divisionColumn + 1) == '-') {
				groupOne.insert(i);
			} else {
				groupTwo.insert(i);
			}
		}
	} 
}

void quickprobs::ColumnRefinement::updateColumnScores(
	const MultiSequence &alignment, std::vector<std::pair<int,float>> &columnScores)
{
	int numSeqs = alignment.count();
	columnScores.resize(alignment.length(), std::pair<int,float>(0,0));

	// get non terminal segment
	std::vector<std::pair<int, int>> nonTerminalFragments(alignment.count());
	for (int i = 0; i < alignment.count(); ++i) { 
		int begin = 0;
		int end = alignment.length() - 1;

		if (config->algorithm.refinement.ignoreTerminalGaps) {
			while (alignment.GetSequence(i)->GetPosition(begin + 1) == '-') {
				++begin;
			}
			while (alignment.GetSequence(i)->GetPosition(end + 1) == '-') {
				--end;
			}	
		}

		nonTerminalFragments[i].first = begin;
		nonTerminalFragments[i].second = end;
	}

	// iterate over columns
	for (int c = 0; c < columnScores.size(); ++c) {
		// iterate over sequences
		columnScores[c].first = c;
		for (int i = 0; i < alignment.count(); ++i) {
			if (c >= nonTerminalFragments[i].first && c <= nonTerminalFragments[i].second) {
				if (alignment.GetSequence(i)->GetPosition(c + 1) == '-') {
					columnScores[c].second += 1.0f;
				}
			}
		}
	}

	std::stable_sort(columnScores.begin(), columnScores.end(), [numSeqs](const std::pair<int,float> &a, const std::pair<int,float>& b){
		return fabs((float)numSeqs / 2 - a.second) > fabs((float)numSeqs / 2 - b.second);
	});

	// erase elements
	auto it = std::remove_if(columnScores.begin(), columnScores.end(), [](const std::pair<int,float> &e){
		return e.second == 0;
	});
	columnScores.erase(it, columnScores.end());
}

