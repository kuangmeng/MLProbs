#pragma once 
#include <string>

#include "Common/ProgramOptions.h"
#include "Alignment/Multiple/Common.h"
#include "Common/Printable.h"
#include "ConfigurationHelper.h"

namespace quickprobs 
{


class Configuration : public Printable
{
public:
	// some I/O configuration

	struct Io {
		std::string input;
		std::string output;
		std::vector<std::string> inputFiles;
		std::vector<std::string> outputFiles;

		std::string annotationFilename;
		bool enableVerbose;
		bool enableAnnotation;
		bool enableClustalWOutput;
		bool enableAlignOrder;
	} io;
	
	// algorithm parameters
	struct Algorithm {
		float posteriorCutoff;
		TreeKind treeKind;
		int smallLimit;
		
		struct Consistency {
			int itertions;
			int smallIterations;
			int largeIterations;
			int iterationsThreshold;
			int numFilterings;
			float saturation;
			
			SelectivityFunction function;
			SelectivityMode mode;
			SelectivityNormalization normalization;
			SelectivityFilter filter;
			float selectivity;
			
			float selfweight;
			float smallSelfweight;
			float largeSelfweight;
			int selfweightThreshold;
			

			SectorCopy copy;

		} consistency;

		struct Refinement {
			RefinementType type;
			int iterations;
			int smallIterations;
			int largeIterations;
			int smallLargeThreshold;
			
			bool acceptanceEntropy;
			bool acceptanceLength;
			float columnFraction;
			int maxDepth;
			bool ignoreTerminalGaps;

			int autosave;
			int seed;
		} refinement;

		float finalSelectivity;
		float finalSaturation;
	
	} algorithm;
	
	quickprobs::AlignmentType type;
	
	struct {
		double gapOpen;
		double gapExtend;
		double temperature;
		std::string matrix;
	} partition;

	struct {
		float gapOpen[2];
		float gapExtend[2];
		std::string name;
		float weight;
	} hmm;

	int stripeCount;
	int stripeLength;
	int binsCount;

	// hardware configuration
	struct {
		int numThreads;
		int refNumThreads;
		int platformNum;
		int deviceNum;
		int64_t memoryLimitMb;
		float gpuMemFactor;
	} hardware;
	
	struct {
		bool usePreciseMath;
		bool useTaylorHmm;
		bool useDoublePartition;
		bool useUnrolledHmm;
		bool divideHmmKernels;
		bool kernelProfiling;
		bool localHmmParams;
		bool localPartitionParams;
	} optimisation;


	Configuration();

	bool parse(int argc, char** argv);

	void setType(quickprobs::AlignmentType type);

	std::string toString();

	std::string printOptions(bool showDeveloperOptions) { return options.toString(showDeveloperOptions); }
	
protected:

	ProgramOptions options;

	void setDefaults();
};

};