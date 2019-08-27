 /*
	QuickProbs version 2.0, October 2015
	Adam Gudys, Sebastian Deorowicz
	Silesian University of Technology, Gliwice, Poland
	adam.gudys@polsl.pl, sebastian.deorowicz@polsl.pl

	GPL version 3.0 applies.
*/
#include <exception>
#include "Common/Log.h"
#include "Common/MemoryTools.h"
#include "Common/rank.h"
#include "Common/deterministic_random.h"

#include "KernelAlignment/Multiple/KernelMSA.h"

using namespace std;

int main(int argc, char* argv[])
{
	std::shared_ptr<quickprobs::KernelMSA> msa;

	TIMER_CREATE(timer);
	
	TIMER_START(timer);
	try {
		TIMER_CREATE(configTimer);
		TIMER_START(configTimer);
		Log::getInstance(Log::LEVEL_NORMAL).enable();
	
		quickprobs::KernelMSA::printWelcome();

		auto config = std::shared_ptr<quickprobs::Configuration>(new quickprobs::Configuration());
		if (argc == 1 || !config->parse(argc, argv)) {
			LOG_NORMAL
				<< "Usage:" << endl	
				<< "\t quickprobs [OPTION]... [infile]..." << endl << endl
				<< config->printOptions(false) << endl
				<< "Following OpenCL devices have been found:" << endl 
				<< clex::OpenCL::listDevices(clex::OpenCL::ANY_DEVICE) << endl << endl; 	
		}
		else {

			if (config->io.enableVerbose) {
				Log::getInstance(Log::LEVEL_DEBUG).enable();
			}

			std::shared_ptr<clex::OpenCL> cl = nullptr;

			if (config->hardware.platformNum >= 0 && config->hardware.deviceNum >= 0) {

				cl = std::shared_ptr<clex::OpenCL>(new clex::OpenCL(
					clex::OpenCL::ANY_DEVICE, config->hardware.platformNum, config->hardware.deviceNum, config->optimisation.kernelProfiling));


				// alter parameters
				if (cl->mainDevice->info->vendor == clex::NVidia) {
					config->optimisation.divideHmmKernels = true;
					config->optimisation.localHmmParams = true;
					config->optimisation.localPartitionParams = true;
				}


				TIMER_STOP(configTimer);

				msa = shared_ptr<quickprobs::KernelMSA>(new quickprobs::KernelMSA(cl, config));
				StatisticsProvider globalStats;

				for (int i = 0; i < config->io.inputFiles.size(); ++i) {
					msa->operator()(config->io.inputFiles[i], config->io.outputFiles[i]);
					globalStats.addStats(*msa);
					msa->reset();
				}

				TIMER_STOP(timer);
			}
		}		
	}
	catch (std::runtime_error &ex) {
		std::cerr << ex.what() << endl;
		return -1;
	}

	return 0;
}
