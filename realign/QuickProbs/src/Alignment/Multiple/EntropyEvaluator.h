#pragma once
#include "AlignmentEvaluatorBase.h"
#include "AminoAcidProperties.h"

namespace quickprobs {

class EntropyEvaluator
{
public:
	virtual float operator()(const MultiSequence& profile);
	virtual float operator()(const MultiSequence& profile, int columnNumber);

protected:
	AminoAcidProperties props;

	static const char* alphabet;
	static const int alphabetSize;
};

}