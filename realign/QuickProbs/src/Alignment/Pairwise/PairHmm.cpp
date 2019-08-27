#include "PairHmm.h"
#include <cassert>

PairHmm::PairHmm(
	int numMatchStates, 
	int numInsertStates, 
	const float* gapOpen, 
	const float* gapExtend, 
	float* initDistrib, 
	float* emitSingle, 
	float* emitPairs, 
	std::string alphabet) : 
	numInsertStates(numInsertStates), numMatchStates(numMatchStates), numStates(numMatchStates + 2 * numInsertStates),
	initDistrib(initDistrib), emitSingle(emitSingle), emitPairs(emitPairs),
	alphabet(alphabet), symbolsCount(alphabet.size())
{
	trans = new float[numStates * numStates];
	
	trans[0] = 1;
	for (int i = 0; i < numInsertStates; ++i) {
		trans[2 * i + 1] = gapOpen[i];
		trans[2 * i + 2] = gapOpen[i];
		trans[0] -= 2 * gapOpen[i]; 
		trans[(2 * i + 1) * numStates + 2 * i + 1] = gapExtend[i];
		trans[(2 * i + 2) * numStates + 2 * i + 2] = gapExtend[i];
		trans[(2 * i + 1) * numStates + 2 * i + 2] = 0;
		trans[(2 * i + 2) * numStates + 2 * i + 1] = 0;
		trans[(2 * i + 1) * numStates] = 1 - gapExtend[i];
		trans[(2 * i + 2) * numStates] = 1 - gapExtend[i];

		assert(trans[0] > 0);
	}
}

