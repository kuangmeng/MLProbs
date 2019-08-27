#include <memory>
#include <algorithm>
#include <omp.h>

#include "Hardware/Kernel.h"
#include "Common/MemoryTools.h"
#include "Common/Log.h"
#include "Common/rank.h"
#include "Common/deterministic_random.h"

#include "Alignment/Multiple/SLinkTree.h"
#include "Alignment/Multiple/ClusterTree.h"
#include "Alignment/Multiple/PhylipTree.h"
#include "Alignment/Multiple/RandomRefinement.h"
#include "Alignment/Multiple/ColumnRefinement.h"
#include "Alignment/Multiple/ScoringRefinement.h"
#include "Alignment/Multiple/TreeRefinement.h"
#include "Alignment/Multiple/ParallelProbabilisticModel.h"

#include "KernelAlignment/KernelFactory.h"
#include "QuickPosteriorStage.h"
#include "QuickConsistencyStage.h"
#include "KernelMSA.h"


#undef VERSION

using namespace std;
using namespace quickprobs;

/// <summary>
/// See declaration for all the details.
/// </summary>
KernelMSA::KernelMSA(std::shared_ptr<clex::OpenCL> cl, std::shared_ptr<Configuration> config)
	: ExtendedMSA(), cl(cl)
{
	TIMER_CREATE(timer);
	TIMER_START(timer);
	
	this->config = config;
	this->version =  std::to_string(QUICKPROBS_VERSION);

	int numCores = omp_get_num_procs();
	LOG_NORMAL << "Detected " << numCores << " CPU cores" << endl;

	if (config->hardware.numThreads <= 0){
		config->hardware.numThreads = numCores;
		LOG_NORMAL << "Automatically set " << config->hardware.numThreads << " OpenMP threads." << endl;
	} else {
		LOG_NORMAL <<"Statically set " << config->hardware.numThreads<<" OpenMP threads." << endl << endl;
	}

	// reassign stages pointer to their kernel variants if OpenCL device specified
	if (cl != nullptr) {
		LOG_NORMAL << "Configuring OpenCL for " << cl->mainDevice->info->deviceName << "." << endl;
		KernelFactory::instance(cl).setFastMath(!config->optimisation.usePreciseMath);
		posteriorStage = std::shared_ptr<PosteriorStage>(new QuickPosteriorStage(cl, config));
		LOG_MEM("After posterior creation");
		consistencyStage = std::shared_ptr<ConsistencyStage>(new QuickConsistencyStage(cl, config));
		LOG_MEM("After consistency creation");
		
	} else {
		LOG_NORMAL << "OpenCL platform and device not specified - using pure CPU version." << endl; 
		posteriorStage = std::shared_ptr<PosteriorStage>(new PosteriorStage(config));
		consistencyStage = std::shared_ptr<ConsistencyStage>(new ConsistencyStage(config));
	}
	
	constructionStage = std::shared_ptr<ConstructionStage>(new ConstructionStage(config));
	
	if (config->algorithm.refinement.type == RefinementType::Column) {
		refinementStage = std::shared_ptr<RefinementBase>(new ColumnRefinement(config, constructionStage));
	} else if (config->algorithm.refinement.type == RefinementType::Tree) {
		refinementStage = std::shared_ptr<RefinementBase>(new TreeRefinement(config, constructionStage));
	} else {
		refinementStage = std::shared_ptr<RefinementBase>(new RandomRefinement(config, constructionStage));
	}

	refinementStage->registerObserver(this);

	TIMER_STOP(timer);
	STATS_WRITE("time.0.1-initialisation", timer.seconds());
}

