#pragma once
#include "PairHmm.h"

class ProteinHmm5 : public PairHmm
{
public:
	ProteinHmm5();

	ProteinHmm5(const float *gapOpen, const float *gapExtend);
};