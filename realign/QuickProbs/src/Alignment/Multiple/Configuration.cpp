#include "Configuration.h"
#include "Alignment/Multiple/Common.h"

#include <iostream>
#include <fstream>

#ifdef WIN32
#include "Common/dirent.h" 
#undef min
#undef max
#else 
#include <dirent.h>
#include <sys/stat.h>
#endif
	
using namespace std;
using namespace quickprobs;


Configuration::Configuration() : options("quickprobs.exe")
{	
	setDefaults();

	// add non-developer parameters

	options.addPositional<std::string>("infile", "input file name");
	options.add<std::string>("outfile,o", "output file name (STDOUT by default)", "", false);
	options.addSwitch("clustalw,l", "use CLUSTALW output format instead of FASTA format", false);
	options.addSwitch("verbose,v", "report progress while aligning", false);
	
	options.add<int>("con-iters,c", "number of consistency repetitions", algorithm.consistency.itertions, false);
	options.add<int>("ref-count,r", "number of iterative refinement passes", algorithm.refinement.iterations, false);
	
	options.addSwitch("nucleotide,n","run QuickProbs in the nucleotide mode", false);

	options.add<int>("num-threads,t", "number of threads (detect automatically if not specified)", 0, false);
	options.add<int>("platform,p", "OpenCL platform id (use CPU mode if not specified)", hardware.platformNum, false);
	options.add<int>("device,d", "OpenCL device id (use CPU mode if not specified)", hardware.deviceNum, false);


	// add developer options
	options.add<float>("cutoff", "", algorithm.posteriorCutoff, true);
	options.add<float>("hmm-weight", "", hmm.weight, true);
	options.add<string>("partition-matrix", "", partition.matrix, true);
	options.add<string>("tree", "guide tree kind", algorithm.treeKind.toString(), true);
	options.add<int>("small-limit", "", algorithm.smallLimit, true);

	options.add<int>("con-postprocess,f", "", algorithm.consistency.numFilterings, true);
	options.add<float>("con-selectivity,x", "", algorithm.consistency.selectivity, true);
	options.add<string>("con-mode", "", algorithm.consistency.mode.toString(), true);
	options.add<string>("con-function", "", algorithm.consistency.function.toString(), true);
	options.add<string>("con-norm", "", algorithm.consistency.normalization.toString(), true);
	options.add<string>("con-filter", "", algorithm.consistency.filter.toString(), true);
	options.add<float>("con-selfweight", "", algorithm.consistency.selfweight, true);
	options.add<string>("sector-copy", "", algorithm.consistency.copy.toString(), true);
	
	options.add<float>("sat3", "sequence weight saturation at consistency", algorithm.consistency.saturation, true);
	options.add<float>("sat4", "sequence weight saturation at final refinement", algorithm.finalSaturation, true);

	options.add<int64_t>("mem-limit", "memory limit", hardware.memoryLimitMb, false);
	
	options.add<string>("ref-type", "random/column/tree", algorithm.refinement.type.toString(), true);
	options.add<string>("ref-acceptance", "acceptance criterion in refinement", "length", true);
	options.add<float>("ref-fraction", "column fraction to be used in refinement", true);
	options.add<int>("ref-depth", "refinement recursion depth", true);
	options.add<int>("ref-autosave", "refinement autosave frequency", true);	
	options.addSwitch("ref-nonterminal", "ignore terminal symbols in refinement", true);
	options.add<int>("ref-seed", "refinement RNG seed", 0, true);
	options.add<int>("ref-threads", "refinement threads", 0, true);

	options.add<bool>("kernel-profiling", "low level kernel profiling", true);
	options.add<bool>("precise-math", "use precise math", true);
	options.add<bool>("double-partition", "use double precision instead of float in partition part", true);
	options.add<bool>("taylor-hmm", "use taylor expansion in HMM part", true);
	options.add<bool>("unrolled-hmm", "use unrolling in HMM part", true);
	options.add<bool>("divided-hmm", "divide HMM kernels", true);
	options.add<bool>("local-hmm-params", "HMM parameters in local memory", true);
	options.add<bool>("local-partition-params", "partition function parameters in local memory", true);
	options.add<float>("gpu-mem-factor", "", 1.0, true);
}

void Configuration::setDefaults()
{
	// some I/O configuration
	io.enableVerbose = false;
	io.enableAnnotation = false;
	io.enableClustalWOutput = false;
	io.enableAlignOrder = false;
	io.annotationFilename = "";
	io.input = "";
	io.output = "";
	
	// algorithm parameters
	hmm.weight = 0.5f;
	partition.matrix = "Vtml200";

	algorithm.treeKind = TreeKind::UPGMA;
	algorithm.posteriorCutoff = 0.01f;

	algorithm.consistency.itertions = -1;
	algorithm.consistency.smallIterations = 2;
	algorithm.consistency.largeIterations = 1;
	algorithm.consistency.iterationsThreshold = 50;

	algorithm.consistency.numFilterings = -1;
	algorithm.consistency.selectivity = 200;
	algorithm.consistency.mode = SelectivityMode::Subtree;
	algorithm.consistency.normalization = SelectivityNormalization::No;
	algorithm.consistency.function = SelectivityFunction::Max;
	algorithm.consistency.filter = SelectivityFilter::Deterministic;

	algorithm.consistency.saturation = 1e-6;

	algorithm.consistency.selfweight = -1.0;
	algorithm.consistency.smallSelfweight = 3.0;
	algorithm.consistency.largeSelfweight = 3.0;
	algorithm.consistency.selfweightThreshold = 200;

	algorithm.consistency.copy = SectorCopy::Full;

	algorithm.finalSelectivity = std::numeric_limits<float>::max();	
	algorithm.finalSaturation = 1e-6;
	
	algorithm.refinement.type = RefinementType::Column;
	algorithm.refinement.iterations = -1;
	algorithm.refinement.smallIterations = 30;
	algorithm.refinement.largeIterations = 200;
	algorithm.refinement.smallLargeThreshold = 200;
	
	algorithm.refinement.columnFraction = 1.0;
	algorithm.refinement.acceptanceLength = true;
	algorithm.refinement.acceptanceEntropy = false;
	algorithm.refinement.maxDepth = 0;
	algorithm.refinement.autosave = std::numeric_limits<int>::max();
	algorithm.refinement.ignoreTerminalGaps = false;

	this->setType(AlignmentType::PROTEIN);
	
	stripeCount = 8;
	stripeLength = 8;
	binsCount = 4;

	// hardware configuration
	hardware.numThreads = 0;
	hardware.refNumThreads = 0;
	hardware.deviceNum = -1;
	hardware.platformNum = -1;
	hardware.memoryLimitMb = static_cast<int64_t>(55e3);
	hardware.gpuMemFactor = 1.0f;

	// optimisation parameters
	optimisation.useDoublePartition = true;
	optimisation.divideHmmKernels = false;
	optimisation.localHmmParams = false;
	optimisation.localPartitionParams = false;
	
	optimisation.usePreciseMath = false;
	optimisation.useTaylorHmm = true;
	optimisation.kernelProfiling = false;
	optimisation.useUnrolledHmm = false;
}

bool Configuration::parse(int argc, char** argv)
{
	try {
		if (!options.parse(argc, argv)) {
			return false;
		}

		// I/O
		options.get("infile", io.input);
		options.get("outfile", io.output);
		options.get("annot", io.annotationFilename);
		options.get("clustalw", io.enableClustalWOutput);
		options.get("verbose", io.enableVerbose);
		options.get("alignment-order", io.enableAlignOrder);

		options.get("cutoff", algorithm.posteriorCutoff);
		options.get("hmm-weight", hmm.weight);
		options.get("partition-matrix", partition.matrix);
		options.get("small-limit", algorithm.smallLimit);

		// algorithm parameters
		string v;
		options.get("ref-type", v);
		algorithm.refinement.type = RefinementType(v);

		options.get("ref-count", algorithm.refinement.iterations);
		options.get("ref-fraction", algorithm.refinement.columnFraction);
		options.get("ref-depth", algorithm.refinement.maxDepth);
		options.get("ref-autosave", algorithm.refinement.autosave);
		options.get("ref-nonterminal", algorithm.refinement.ignoreTerminalGaps);
		options.get("ref-seed", algorithm.refinement.seed);
	
		if (io.enableVerbose) {
			algorithm.refinement.autosave = 5;
		}

		string acceptance;
		options.get("ref-acceptance", acceptance);
		if (acceptance == "length") {
			algorithm.refinement.acceptanceLength = true;
			algorithm.refinement.acceptanceEntropy = false;
		} else if (acceptance == "entropy") {
			algorithm.refinement.acceptanceLength = false;
			algorithm.refinement.acceptanceEntropy = true;
		} else if (acceptance == "length_entropy") {
			algorithm.refinement.acceptanceLength = true;
			algorithm.refinement.acceptanceEntropy = true;
		} else {
			algorithm.refinement.acceptanceLength = false;
			algorithm.refinement.acceptanceEntropy = false;
		}

		options.get("con-iters", algorithm.consistency.itertions);
		options.get("con-postprocess", algorithm.consistency.numFilterings);
		options.get("con-selectivity", algorithm.consistency.selectivity);
		options.get("con-selfweight", algorithm.consistency.selfweight);

		string t;
		options.get("con-norm", t);
		algorithm.consistency.normalization = SelectivityNormalization(t);
		options.get("con-mode", t);
		algorithm.consistency.mode = SelectivityMode(t);
		options.get("con-function", t);
		algorithm.consistency.function = SelectivityFunction(t);
		options.get("con-filter", t);
		algorithm.consistency.filter = SelectivityFilter(t);
		options.get("sector-copy", t);
		algorithm.consistency.copy = SectorCopy(t);
		options.get("tree", t);
		algorithm.treeKind = TreeKind(t);


		options.get("sat3", algorithm.consistency.saturation);
		options.get("sat4", algorithm.finalSaturation);
	
		bool nucleotideMode;
		options.get("nucleotide", nucleotideMode);
		this->setType(nucleotideMode ? AlignmentType::NUCLEOTIDE : AlignmentType::PROTEIN);

		// hardware parameters
		options.get("num-threads", hardware.numThreads);
		options.get("platform", hardware.platformNum);
		options.get("device", hardware.deviceNum);
		options.get("mem-limit", hardware.memoryLimitMb);

		DIR *inputDir = opendir(io.input.c_str());
		DIR *outputDir = opendir(io.output.c_str());
		struct dirent *ent;
		if ((inputDir != NULL) && (outputDir != NULL)) {
			while ((ent = readdir(inputDir)) != NULL) {
				struct stat info;
				string inputPath = io.input + "/" + string(ent->d_name);
				stat(inputPath.c_str(), &info);
				// only for regular files
				if (S_ISREG(info.st_mode)) {
					io.inputFiles.push_back(inputPath);
					io.outputFiles.push_back(io.output + "/" + string(ent->d_name));
				}
			}
			closedir(inputDir);
			closedir(outputDir);
		}
		else {
			io.inputFiles.push_back(io.input);
			io.outputFiles.push_back(io.output);
		}

		/* std::filesystem variant
		fs::path inputPath(io.input);
		fs::path outputPath(io.output);
		if (fs::is_directory(inputPath)) {
			if (fs::is_directory(outputPath)) {
				for (auto it = fs::directory_iterator(inputPath); it != fs::directory_iterator(); ++it) {
					auto file = it->path().filename();
					io.inputFiles.push_back(it->path().string());
					io.outputFiles.push_back(outputPath.string() + "//" + file.c_str());
				}
			}
		} else {
			io.inputFiles.push_back(io.input);
			io.outputFiles.push_back(io.output);
		}
		*/

		// optimisation parameters
		//options.get("kernel-profiling", optimisation.kernelProfiling);
		options.get("precise-math", optimisation.usePreciseMath);
		options.get("taylor-hmm", optimisation.useTaylorHmm);
		//options.get("unrolled-hmm", optimisation.useUnrolledHmm);
		options.get("divided-hmm", optimisation.divideHmmKernels);
	//	options.get("double-partition", optimisation.useDoublePartition);
		options.get("local-hmm-params", optimisation.localHmmParams);
		options.get("local-partition-params", optimisation.localPartitionParams);
		options.get("ref-threads", hardware.refNumThreads);
		options.get("gpu-mem-factor", hardware.gpuMemFactor);
	}
	catch (std::runtime_error& error) {
		return false;
	}
	
	return true;
}

void quickprobs::Configuration::setType(quickprobs::AlignmentType type)
{
	this->type = type;
	
	if (type == quickprobs::AlignmentType::PROTEIN) {
		hmm.name = "ProteinHmm5";
		hmm.gapExtend[0] = 0.3965826333f;
		hmm.gapExtend[1] = 0.8988758326f;
		hmm.gapOpen[0] = 0.0119511066f;
		hmm.gapOpen[1] = 0.008008334786f;
		
		if (partition.matrix == "MiqsFP") {
			partition.gapExtend = -1.31272;
			partition.gapOpen = -22.9675;
			partition.temperature = 5.02492;
		} else if (partition.matrix == "Vtml200") {
	//		partition.gapExtend = -1.5;
	//		partition.gapOpen = -22.15;
	//		partition.temperature = 4.95;

	//		partition.gapExtend = -1;
	//		partition.gapOpen = -22;
	//		partition.temperature = 5;

			partition.gapExtend = -1.30113;
			partition.gapOpen = -25.3549;
			partition.temperature = 5.6007;
		} else if (partition.matrix == "Gonnet160") {
			partition.gapExtend = -1;
			partition.gapOpen = -22;
			partition.temperature = 5;
		}
		
		// optimised parameters
	/*	hmm.gapExtend[0] = 0.499997f;
		hmm.gapExtend[1] = 0.893421f;
		hmm.gapOpen[0] = 0.0587073f;
		hmm.gapOpen[1] = 0.00852711f;
		partition.gapExtend = -1.79459f;
		partition.gapOpen = -27.5605f;
		partition.temperature = 6.0269f;
*/
	} else {
		partition.gapOpen = -414;
		partition.gapExtend = -27;
		partition.temperature = 100.0;
		partition.matrix = "Hoxd70";

		hmm.gapOpen[0] = 0.0119511066f;
		hmm.gapOpen[1] = 0.008008334786f;
		hmm.gapExtend[0] = 0.3965826333f;
		hmm.gapExtend[1] = 0.8988758326f;
		hmm.name = "NucleotideHmmGTR5";
		//hmm.name = "ProteinHmm5";
	}
}

std::string quickprobs::Configuration::toString()
{
	std::ostringstream ss;
	ss << "I/O:" << endl
		<< "inputFilename=" << io.input << endl
		<< "outputFilename=" << io.output << endl
		<< "annotationFilename=" << io.annotationFilename << endl
		<< "enableClustalWOutput=" << io.enableClustalWOutput << endl
		<< "enableVerbose=" << io.enableVerbose << endl
		<< "enableAlignOrder=" << io.enableAlignOrder << endl;

	ss << endl << "HARDWARE:" << endl
		<< "memoryLimitMb=" << hardware.memoryLimitMb << endl
		<< "gpuMemFactor=" << hardware.gpuMemFactor << endl; 

	ss << endl << "Algorithm:" << endl
		<< "partition.matrix=" << partition.matrix << endl
		<< "partition.gapOpen=" << partition.gapOpen << endl
		<< "partition.gapExtend=" << partition.gapExtend << endl
		<< "partition.temperature=" << partition.temperature << endl
		<< "hmm.weight=" << hmm.weight << endl
		<< "degenerateTree=" << algorithm.treeKind.toString() << endl
		<< "consistency.itertions=" << algorithm.consistency.itertions << endl
		<< "consistency.numFilterings=" << algorithm.consistency.numFilterings << endl
		<< "consistency.selectivity=" << algorithm.consistency.selectivity << endl
		<< "consistency.saturation=" << algorithm.consistency.saturation << endl	
		<< "consistency.normalization=" << algorithm.consistency.normalization.toString() << endl
		<< "consistency.mode=" << algorithm.consistency.mode.toString() << endl
		<< "consistency.function=" << algorithm.consistency.function.toString() << endl
		<< "consistency.filter=" << algorithm.consistency.filter.toString() << endl
		<< "consistency.selfweight=" << algorithm.consistency.selfweight << endl
		<< "consistency.sectorCopy=" << algorithm.consistency.copy.toString() << endl
		<< "finalSelectivity=" << algorithm.finalSelectivity << endl
		<< "finalSaturation=" << algorithm.finalSaturation << endl
		<< "posteriorCutoff=" << algorithm.posteriorCutoff << endl
		<< "smallSeqlimit=" << algorithm.smallLimit << endl
		<< "refinement.type=" << algorithm.refinement.type.toString() << endl
		<< "refinement.iterations=" << algorithm.refinement.iterations << endl
		<< "refinement.columnFraction=" << algorithm.refinement.columnFraction << endl
		<< "refinement.maxDepth=" << algorithm.refinement.maxDepth << endl
		<< "refinement.acceptanceLength=" << algorithm.refinement.acceptanceLength << endl
		<< "refinement.acceptanceEntropy=" << algorithm.refinement.acceptanceEntropy << endl;;

	ss << endl << "Optimisation:" << endl 
	//	<< "kernelProfiling=" << optimisation.kernelProfiling << endl
		<< "usePreciseMath=" << optimisation.usePreciseMath << endl
		<< "useTaylorHmm=" << optimisation.useTaylorHmm << endl
		<< "divideHmmKernels=" << optimisation.divideHmmKernels << endl
		<< "useDoublePartition=" << optimisation.useDoublePartition << endl
		<< "localHmmParams=" << optimisation.localHmmParams << endl
		<< "localPartitionParams=" << optimisation.localPartitionParams << endl;
	

	return ss.str();
}
