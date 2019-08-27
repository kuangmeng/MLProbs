#pragma once
#include "PairHmm.h"

class PfamHmm5 : public PairHmm
{
public:
	PfamHmm5();

	PfamHmm5(const float *gapOpen, const float *gapExtend);
};