#include "NucleotideHmmGTR5.h"

namespace {
	float initDistribDefault[] = { 0.6814756989f, 8.615339902e-05f, 8.615339902e-05f, 0.1591759622f, 0.1591759622f };
	float gapOpenDefault[] = { 0.0119511066f, 0.0119511066f, 0.008008334786f, 0.008008334786f };
	float gapExtendDefault[] = { 0.3965826333f, 0.3965826333f, 0.8988758326f, 0.8988758326f };

	float emitSingleDefault[] = {
		0.35, 0.15, 0.25, 0.25,
	};

	float emitPairsDefault[] = {
		0.1109,    0.0324,    0.0661,    0.0545,
		0.0744,    0.0504,    0.0450,    0.0741,
		0.0947,    0.0289,    0.0760,    0.0489,
		0.0725,    0.0438,    0.0429,   0.0844
	};
}

NucleotideHmmGTR5::NucleotideHmmGTR5() 
	: PairHmm(1, 2, gapOpenDefault, gapExtendDefault, initDistribDefault, emitSingleDefault, emitPairsDefault, "ACGT") 
{
}
