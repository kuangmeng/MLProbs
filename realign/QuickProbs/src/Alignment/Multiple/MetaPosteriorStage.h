#pragma once
#include "QuickPosteriorStage.h"

namespace quickprobs
{

class MetaPosteriorStage : public QuickPosteriorStage
{
public:
	MetaPosteriorStage(
		std::shared_ptr<clex::OpenCL> cl, 
		std::shared_ptr<Configuration> config);

protected:

	virtual void run(
		MultiSequence& sequences, 
		std::vector<std::vector<float>>& distances,
		std::vector<std::vector<SparseMatrixType*>>& matrices);

};

}