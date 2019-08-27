#pragma once

#include "IPairwiseAligner.h"

class NeedlemanWunsch : public IPairwiseAligner
{
public:
	const int *substitution;
	int symbolsCount;
	int gi;
	int ge;

	
	NeedlemanWunsch(const int *substitution, int symbolsCount, int gi, int ge) :
		substitution(substitution), symbolsCount(symbolsCount), gi(gi), ge(ge) {}

	virtual int operator()(
		const char *seq1, 
		const char *seq2,
		int seq1Length,
		int seq2Length);

	virtual int operator()(
		const char *seq1, 
		const char *seq2,
		int seq1Length,
		int seq2Length,
		std::vector<int>& scores,
		std::vector<int>& backtrack);
};