#include <algorithm>

#include "ContiguousMultiSequence.h"
#include "Sequence.h"
#include "MultiSequence.h"
#undef max
#undef min

using namespace quickprobs;

ContiguousMultiSequence::ContiguousMultiSequence(const MultiSequence& reference)
{
	this->sequences = reference.sequences;
	int numSeq = reference.count();
	int totalSize = 0;
	maxLength = 0;
	minLength = std::numeric_limits<unsigned int>::max();

	// iterate for the first time just to compute size of data
	for (int i = 0; i < numSeq; ++i) {
		totalSize += reference.GetSequence(i)->GetLength() + 1; 
	}
	totalSize++; // for ending zero

	// iterate second time
	lengths.resize(numSeq);
	offsets.resize(numSeq);
	sortingMap.resize(numSeq);
	sequenceData.resize(totalSize);
	
	auto dataPtr = sequenceData.begin();
	int offset = 0;

	for (int i = 0; i < numSeq; ++i) {
		auto ref = reference.GetSequence(i);
		auto seq = encode(*ref, nullptr);
		unsigned int length = seq->GetLength();

		std::copy(seq->getIterator(), seq->getIterator() + length + 1, dataPtr);

		maxLength = std::max(maxLength, length);
		minLength = std::min(minLength, length);

		lengths[i] = length;
		offsets[i] = offset;
		offset += length + 1;
		dataPtr += length + 1; 	
	}

	*dataPtr = 0; // add finishing zero
	
	int i = 0;
	std::generate(sortingMap.begin(), sortingMap.end(), [&i]()->unsigned int {
		return i++;
	});

	std::stable_sort(sortingMap.begin(), sortingMap.end(), [this](unsigned int a, unsigned b)->bool {
		return lengths[a] < lengths[b];
	});
}

std::unique_ptr<Sequence> ContiguousMultiSequence::encode(const Sequence& src, const int* map)
{
	auto out = std::unique_ptr<Sequence>(src.Clone());
	
	const char* srcData = src.getData();
	char* outData = out->getData();

	for (int i = 1; i <= src.GetLength(); i++) // omit @ symbol at the beginning of the sequence
	{
		outData[i] = srcData[i] - 'A';
	}

	return out;
}


