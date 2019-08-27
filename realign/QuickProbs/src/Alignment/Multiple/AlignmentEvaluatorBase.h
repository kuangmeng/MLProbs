#pragma once
#include "DataStructures/MultiSequence.h"
#include <set>
#include <algorithm>

namespace quickprobs {

class AlignmentEvaluatorBase 
{
public:
	virtual float operator()(const MultiSequence& alignment) = 0;
};

}