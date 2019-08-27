#pragma once
#include "Alignment/DataStructures/MultiSequence.h"

namespace quickprobs
{

class IRefinementObserver 
{
public:
	virtual void iterationDone(const MultiSequence& alignment, int iteration) = 0;
};

}