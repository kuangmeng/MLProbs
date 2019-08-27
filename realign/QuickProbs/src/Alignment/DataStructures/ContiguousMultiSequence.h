#pragma once
#include <vector>
#include <memory>
#include "ISequenceSet.h"

namespace quickprobs 
{
class Sequence;
class MultiSequence;

class ContiguousMultiSequence : public ISequenceSet
{
public:
	std::vector<char> sequenceData;
	std::vector<unsigned int> offsets;
	std::vector<unsigned int> lengths;
	std::vector<unsigned int> sortingMap;
	std::vector<quickprobs::Sequence*> sequences;

	unsigned int maxLength;
	unsigned int minLength;

	virtual ::size_t count() const { return offsets.size(); }

	/// Retrieves a sequence from the MultiSequence object.
	virtual Sequence* GetSequence(int i){ return sequences[i]; }

	/// Retrieves a sequence from the MultiSequence object (const version).
	virtual const Sequence* GetSequence(int i) const { return sequences[i]; }

	ContiguousMultiSequence(const MultiSequence& reference);

	static std::unique_ptr<Sequence> encode(const Sequence& src, const int* map);
};

};