#include "EntropyEvaluator.h"
#include "AminoAcidProperties.h"
#include "Common/mathex.h"

#include <set>
#include <cstring>
#include <bitset>

using namespace std;
using namespace quickprobs;

const char* EntropyEvaluator::alphabet = "ARNDCQEGHILKMFPSTWYV";
const int EntropyEvaluator::alphabetSize = 20;

float quickprobs::EntropyEvaluator::operator()(const MultiSequence& profile)
{
	float score = 0;
		
	// iterate over alignment columns
	for (int c = 0; c < profile.length(); ++c) {
		score += (*this)(profile, c + 1);
	}

	return score;
}

float quickprobs::EntropyEvaluator::operator()(const MultiSequence& profile, int columnNumber)
{
	unsigned int commonProperties = 0xffffffff;
	unsigned int allProperties = 0x0;
	float lambda = 1.0f / std::log2(alphabetSize);
	
	// for each symbol in column
	std::vector<float> histogram(256, 0.5); // initialise with 0.5
	float entropy = 0;
	int gapCount = 0;
	int symbolsCount = alphabetSize / 2;

	// get histogram
	for (auto i = 0; i != profile.count(); ++i) {
		char s = profile.GetSequence(i)->GetPosition(columnNumber);
		if (s != '-') {
			histogram[s] += 1.0f; // add 1 for each occurrence
			++symbolsCount; 

			commonProperties &= props.get(s);
			allProperties |= props.get(s);
		} else {
			++gapCount;
		}
	}

	// property score
	std::bitset<10> commonBits(commonProperties);
	std::bitset<10> allBits(allProperties);
	float propScore = (float)commonBits.count() + 10.0f - (float)allBits.count();
	propScore /= commonBits.size();

	// entropy score
	for (int i = 0; i < alphabetSize; ++i) {
		char s = alphabet[i];
		double ps = (double)histogram[s] / symbolsCount;
		entropy -= lambda * (float)(ps * std::log2(ps));
	}

	// gap score
	float gapScore = (float)gapCount / profile.count();

	// total score
	float score = (1 - entropy) * propScore * (1 - gapScore);
	return score;
}

