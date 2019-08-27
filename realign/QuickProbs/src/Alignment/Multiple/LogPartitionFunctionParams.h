#pragma once
#include "IPartitionFunctionParams.h"
#include "Common.h"
#include "ScoreType.h"
#include "Alignment/Pairwise/Scoring.h"

#include <CL/cl.h>
#include <cstring>


namespace quickprobs
{

class LogPartitionFunctionParams : public IPartitionFunctionParams
{
public:
	struct Raw {
		float termGapOpen; 
		float termGapExtend; 
		float gapOpen; 
		float gapExt;
		float subMatrix[26 * 26];
	} raw;

	inline LogPartitionFunctionParams(const Scoring<double>& scoring, double temperature);
	::size_t sizeInBytes() const { return sizeof(raw); }
	void* data() const { return (void*)(&raw); }
};

LogPartitionFunctionParams::LogPartitionFunctionParams(const Scoring<double>& scoring, double temperature)
{
	double beta = 1.0 / temperature;
	std::memset(raw.subMatrix, LOG_ZERO, sizeof(raw.subMatrix));

	for (int i = 0; i < scoring.symbolsCount - 1; i++) { // ignore last element
		int out_i = scoring.alphabet[i] - 'A';
		for (int j = 0; j <= i; j++) {
			int out_j = scoring.alphabet[j] - 'A';

			double value = beta * scoring.get(i, j);
			raw.subMatrix[out_i * 26 + out_j] = raw.subMatrix[out_j * 26 + out_i] = (float)value;
		}
	}

	raw.gapOpen = (float)(beta * scoring.gi);
	raw.gapExt = (float)(beta * scoring.ge);
	raw.termGapOpen = (float)(beta * 0);
	raw.termGapExtend = (float)(beta * 0);
}


};