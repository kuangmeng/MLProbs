#pragma once
#include "PairHmm.h"
#include "ProteinHmm3.h"
#include "ProteinHmm5.h"
#include "PfamHmm5.h"
#include "NucleotideHmmGTR5.h"
#include <memory>
#include <stdexcept>

class PairHmmFactory {
public:
	static std::unique_ptr<PairHmm> create(std::string name, const float *gapOpen, const float* gapExtend) {

		PairHmm* out = nullptr;

		if		(name == "ProteinHmm3")			{ out = new ProteinHmm3(); }
		else if (name == "ProteinHmm5")			{ out = new ProteinHmm5(gapOpen, gapExtend); }
		else if (name == "PfamHmm5")			{ out = new PfamHmm5(gapOpen, gapExtend);	 }
		else if (name == "NucleotideHmmGTR5")	{ out = new NucleotideHmmGTR5(); }
		else throw 
			std::runtime_error("PairHmmFactory::create(): Unknown HMM type.");

		return std::unique_ptr<PairHmm>(out);
		
	}
};