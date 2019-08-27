#include "KernelSet.h"

#include "KernelAlignment/KernelRepository.h"
#include "KernelAlignment/KernelFactory.h"

using namespace quickprobs;

void KernelSet::fillAll(const Configuration& config, std::shared_ptr<clex::OpenCL> cl)
{
	std::vector<std::string> defines;
	std::vector<int> commonFiles;
	std::vector<int> files;
	KernelFactory& factory = KernelFactory::instance(cl);

	if (config.optimisation.useTaylorHmm)			{ defines.push_back("USE_TAYLOR_HMM"); }
	if (!config.optimisation.useDoublePartition)	{ defines.push_back("USE_FLOAT_PARTITION"); }
	if (config.optimisation.divideHmmKernels)		{ defines.push_back("RETURN_GLOBAL"); }
	if (config.optimisation.localHmmParams)			{ defines.push_back("COPY_PROBABILISTIC_PARAMS"); }
	if (config.optimisation.localPartitionParams)	{ defines.push_back("COPY_PARTITION_PARAMS"); }

	defines.push_back("POSTERIOR_CUTOFF=" + std::to_string(config.algorithm.posteriorCutoff));
	defines.push_back("HMM_WEIGHT=" + std::to_string(config.hmm.weight));

	commonFiles.push_back(KernelRepository::Utils);
	commonFiles.push_back(KernelRepository::MemoryTransfers);
	commonFiles.push_back(KernelRepository::ScoreType);
	commonFiles.push_back(KernelRepository::AuxiliaryStructures);
	commonFiles.push_back(KernelRepository::LogarithmUtils);

	// compile classic kernels
	files = commonFiles;
	if (config.optimisation.useDoublePartition) {
		files.push_back(KernelRepository::Partition);
	} else {
		files.push_back(KernelRepository::PartitionLogarithm);
	}

	files.push_back(KernelRepository::Probabilistic);
	files.push_back(KernelRepository::Finalization);

	files.push_back(KernelRepository::MultiplePartition);
	defines.push_back("LOG_UNDERFLOW=50.0f");
	insert(PARTITION_FORWARD, factory.create(files, "MultiplePartition_ComputeForward", defines), 
		config.optimisation.useDoublePartition ? 6 : 3);
	insert(PARTITION_REVERSE, factory.create(files, "MultiplePartition_ComputeReverse", defines), 
		config.optimisation.useDoublePartition ? 8 : 4);
	defines.pop_back();
	files.pop_back();
	

	files.push_back(KernelRepository::MultipleProbabilistic);
	defines.push_back("LOG_UNDERFLOW=7.5f");
	if (config.optimisation.divideHmmKernels) {
		insert(HMM_FORWARD, factory.create(files, "MultipleProbabilistic_ComputeForward", defines), 3);
		insert(HMM_BACKWARD, factory.create(files, "MultipleProbabilistic_ComputeBackward", defines), 2);
		insert(HMM_COMBINE, factory.create(files, "MultipleProbabilistic_Combine", defines), 0);
	} else {
		insert(HMM_ALL, factory.create(files, "MultipleProbabilistic_ComputeAll", defines), 3);
	}
	files.pop_back();

	files.push_back(KernelRepository::MultipleFinalization);
	insert(FINALIZATION, factory.create(files, "MultipleFinalization_Combine", defines), 1);

	files = commonFiles;
	files.push_back(KernelRepository::SparseMatrix);
	files.push_back(KernelRepository::MultipleSparse);
	insert(SPARSE, factory.create(files, "MultipleSparse_Compute", defines), 0);
}

void quickprobs::KernelSet::fillAllLong(const Configuration& config, std::shared_ptr<clex::OpenCL> cl)
{
	std::vector<std::string> defines;
	std::vector<int> commonFiles;
	std::vector<int> files;
	KernelFactory& factory = KernelFactory::instance(cl);
	
	//if (config.optimisation.useTaylorHmm)			{ defines.push_back("USE_TAYLOR_HMM"); }
	if (!config.optimisation.useDoublePartition)		{ defines.push_back("USE_FLOAT_PARTITION"); }
	if (config.optimisation.divideHmmKernels)		{ defines.push_back("RETURN_GLOBAL"); }
	if (config.optimisation.localHmmParams)			{ defines.push_back("COPY_PROBABILISTIC_PARAMS"); }
	if (config.optimisation.localPartitionParams)	{ defines.push_back("COPY_PARTITION_PARAMS"); }

	defines.push_back("POSTERIOR_CUTOFF=" + std::to_string(config.algorithm.posteriorCutoff));
	defines.push_back("HMM_WEIGHT=" + std::to_string(config.hmm.weight));

	commonFiles.push_back(KernelRepository::Utils);
	commonFiles.push_back(KernelRepository::MemoryTransfers);
	commonFiles.push_back(KernelRepository::ScoreType);
	commonFiles.push_back(KernelRepository::AuxiliaryStructures);
	commonFiles.push_back(KernelRepository::LogarithmUtils);
	
	// compile long kernels
	defines.push_back("LONG_KERNELS");
	files = commonFiles;
	if (config.optimisation.useDoublePartition) {
		files.push_back(KernelRepository::Partition_long);
	} else {
		files.push_back(KernelRepository::PartitionLogarithm_long);
	}

	files.push_back(KernelRepository::Probabilistic_long);
	files.push_back(KernelRepository::Finalization_long);

	files.push_back(KernelRepository::MultiplePartition);
	defines.push_back("LOG_UNDERFLOW=50.0f");
	insert(PARTITION_FORWARD, factory.create(files, "MultiplePartition_ComputeForward_long", defines), 
		config.optimisation.useDoublePartition ?  6 : 3);
	insert(PARTITION_REVERSE, factory.create(files, "MultiplePartition_ComputeReverse_long", defines),
		config.optimisation.useDoublePartition ? 8 : 4);
	defines.pop_back();
	files.pop_back();

	files.push_back(KernelRepository::MultipleProbabilistic);
	defines.push_back("LOG_UNDERFLOW=7.5f");
	if (config.optimisation.divideHmmKernels) {
		insert(HMM_FORWARD, factory.create(files, "MultipleProbabilistic_ComputeForward_long", defines), 3);
		insert(HMM_BACKWARD, factory.create(files, "MultipleProbabilistic_ComputeBackward_long", defines), 2);
		insert(HMM_COMBINE, factory.create(files, "MultipleProbabilistic_Combine_long", defines), 0);
	} else {
		insert(HMM_ALL, factory.create(files, "MultipleProbabilistic_ComputeAll_long", defines), 3);
	}
	files.pop_back();

	files.push_back(KernelRepository::MultipleFinalization);
	insert(FINALIZATION, factory.create(files, "MultipleFinalization_Combine_long", defines), 1);

	files = commonFiles;
	files.push_back(KernelRepository::SparseMatrix_long);
	files.push_back(KernelRepository::MultipleSparse);
	insert(SPARSE, factory.create(files, "MultipleSparse_Compute_long", defines), 0);
}



