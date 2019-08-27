#pragma once
#include <string>
#include <vector>
#include <CL/cl.h>

#include "Alignment/Pairwise/PairHmm.h"

namespace quickprobs
{

struct ProbabilisticParams
{
	cl_float initial[5];			// holds the initial probabilities for each state
	cl_float trans[5 * 5];			// holds all state-to-state transition probabilities
	cl_float match[26 * 26];
	cl_float insert[26 * 5];

	void init (const PairHmm& hmm);
	::size_t sizeInBytes() { return sizeof(ProbabilisticParams); }
	void* data() { return (void*)(this); }
};

};