#pragma once
#include <memory>

#include "Hardware/OpenCl.h"
#include "Hardware/Kernel.h"
#include "Alignment/Multiple/ConsistencyStage.h"
#include "Alignment/Multiple/Configuration.h"

#include "RelaxationSector.h"
#include "RelaxationKernelWrapper.h"

#include "Common/Array.h"

namespace quickprobs
{

class QuickConsistencyStage : public ConsistencyStage
{
public:
	QuickConsistencyStage(
		std::shared_ptr<clex::OpenCL> cl,
		std::shared_ptr<Configuration> config);

protected:
	
	std::shared_ptr<clex::OpenCL> cl; 

	std::vector<std::shared_ptr<RelaxationKernelWrapper>> relaxKernels;

	std::shared_ptr<clex::Kernel> kernel_postprocessUnfiltered;
	std::shared_ptr<clex::Kernel> kernel_postprocessFiltered;

	virtual void run(
		const float* seqsWeights,
		ISequenceSet& set, 
		Array<float>& distances,
		Array<SparseMatrixType*>& sparseMatrices);

	virtual void doRelaxationGpu(
		const float* seqsWeights, 
		ContiguousMultiSequence& sequences, 
		Array<float>& distances,
		Array<SparseMatrixType*>& sparseMatrices,
		bool filter);

	std::vector<std::shared_ptr<RelaxationSector>> generateSectors(
		const Array<SparseMatrixType*>& sparseMatrices,
		const ContiguousMultiSequence& sequences,
		const Array<float>& distances,
		std::vector<unsigned int>& sparseOffsets,
		::size_t& maxSectorBytes);

	void printSelectivityHistogram(std::vector<std::shared_ptr<RelaxationSector>> sectors);

	void printDistanceHistogram(Array<float>& distances);

	void printSparseRowsHistogram(Array<SparseMatrixType*>& sparseMatrices);

};

};