#include <iostream>
#include <cassert>
#include <algorithm>

#include "ProbabilisticParams.h"
#include "Pairwise/PairHmm.h"
#include "ScoreType.h"

using namespace std;
using namespace quickprobs;

void quickprobs::ProbabilisticParams::init(const PairHmm& hmm)
{
	
	// create initial and transition probability matrices
	for (int i = 0; i < hmm.numStates; ++i) {
		initial[i] = LOG(hmm.initDistrib[i]);
		for (int j = 0; j < hmm.numStates; ++j)
			trans[i * hmm.numStates + j] = LOG(hmm.trans[i * hmm.numStates + j]);
	}

	// initialise arrays
	float log_e5 = LOG(1e-5);
	float log_e10 = LOG(1e-10);
	std::fill(insert, insert + 26 * 5, log_e5);
	std::fill(match, match + 26 * 26, log_e10);

	// create insertion and match probability matrices
	for (int i = 0; i < hmm.symbolsCount; ++i) {
		int out_i = hmm.alphabet[i] - 'A';

		for (int d = 0; d < hmm.numStates; ++d) {
			insert[out_i * hmm.numStates + d] = LOG(hmm.emitSingle[i]);
		}

		for (int j = 0; j <= i; ++j) {
			int out_j = hmm.alphabet[j] - 'A';
			match[out_i * 26 + out_j] = match[out_j * 26 + out_i] = 
				LOG(hmm.emitPairs[i * hmm.symbolsCount + j]);
		}
	}
}
