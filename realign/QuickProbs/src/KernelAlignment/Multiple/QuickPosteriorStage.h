#pragma once
#include "Hardware/OpenCl.h"
#include "Hardware/Kernel.h"
#include "Alignment/Multiple/PosteriorStage.h"
#include "PosteriorTasksWave.h"
#include "KernelSet.h"

#include <atomic>

namespace quickprobs
{

class QuickPosteriorStage : public PosteriorStage
{
public:
	QuickPosteriorStage(std::shared_ptr<clex::OpenCL> cl, std::shared_ptr<Configuration> config);

protected:
	std::shared_ptr<clex::OpenCL> cl; 

	std::vector<KernelSet> kernelSets;
	
	std::vector<int> maxWorkgroupSizes;

	int64_t sparseBytes;

	int64_t denseBytes;

	double timeCpuProcessing;

	virtual void run(
		ISequenceSet& set,  
		Array<float>& distances,
		Array<SparseMatrixType*>& matrices);

	virtual void computeWaveCpu(
		const PosteriorTasksWave& wave,
		const ContiguousMultiSequence& sequences,
		Array<float>& distances,
		Array<SparseMatrixType*>& matrices,
		std::atomic< ::size_t>& sparseBytes);
};

};