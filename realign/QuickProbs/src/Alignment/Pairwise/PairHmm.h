#pragma once
#include <string>

class PairHmm
{
public:
	const int numStates;
	const int numInsertStates;
	const int numMatchStates;
	
	float* trans;
	const float* initDistrib;
	const float* emitSingle;
	const float* emitPairs;
	
	const std::string alphabet;
	const int symbolsCount;

	PairHmm(
		int numMatchStates, 
		int numInsertStates, 
		const float* gapOpen, 
		const float* gapExtend, 
		float* initDistrib, 
		float* emitSingle, 
		float* emitPairs,
		std::string alphabet);

	~PairHmm() { delete [] trans; }
};
