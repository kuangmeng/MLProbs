#pragma once
#include "IPartitionFunctionParams.h"
#include "Common.h"
#include "Alignment/Pairwise/Scoring.h"
#include "Common/mathex.h"

#include <cstring>

namespace quickprobs
{

template <class partition_t>
class ExpPartitionFunctionParams : public IPartitionFunctionParams
{
public:
	struct Raw {
		partition_t termGapOpen; 
		partition_t termGapExtend; 
		partition_t gapOpen; 
		partition_t gapExt;
		partition_t subMatrix[26 * 26];
	} raw;

	ExpPartitionFunctionParams(const Scoring<double>& scoring, double temperature);
	::size_t sizeInBytes() const { return sizeof(raw); }
	void* data() const { return (void*)(&raw); }
};

template <class partition_t>
ExpPartitionFunctionParams<partition_t>::ExpPartitionFunctionParams(const Scoring<double>& scoring, double temperature)
{
	double beta = 1.0 / temperature;
	std::memset(raw.subMatrix, 0, sizeof(raw.subMatrix));

	for (int i = 0; i < scoring.symbolsCount - 1; i++) { // ignore last element
		int out_i = scoring.alphabet[i] - 'A';
		for (int j = 0; j <= i; j++) {
			int out_j = scoring.alphabet[j] - 'A';

			double value = exp(beta * scoring.get(i, j));
			raw.subMatrix[out_i * 26 + out_j] = raw.subMatrix[out_j * 26 + out_i] = mathex::convertUp<double, partition_t>(value);
		}
	}

	raw.gapOpen = mathex::convertUp<double, partition_t>(exp(beta * scoring.gi));
	raw.gapExt = mathex::convertUp<double, partition_t>(exp(beta * scoring.ge));
	raw.termGapOpen = mathex::convertUp<double, partition_t>(exp(beta * 0));
	raw.termGapExtend = mathex::convertUp<double, partition_t>(exp(beta * 0));
}


};