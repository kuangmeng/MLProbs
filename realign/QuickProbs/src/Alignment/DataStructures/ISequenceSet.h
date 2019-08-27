#pragma once

namespace quickprobs 
{

// forward declaration
class Sequence;

// class definition
class ISequenceSet
{
public:
	virtual Sequence* GetSequence (int i) = 0;
	virtual const Sequence* GetSequence (int i) const = 0;
	virtual ::size_t count() const = 0;
	
	virtual ~ISequenceSet() {}
};

};