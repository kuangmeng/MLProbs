#include "Common/mathex.h"
#include "Common/Log.h"
#include "Common/MemoryTools.h"
#include "Common/deterministic_random.h"

#include "Alignment/Multiple/MemoryPool.h"
#include "Alignment/Multiple/GuideTree.h"
#include "Alignment/DataStructures/ContiguousMultiSequence.h"
#include "Alignment/DataStructures/SparseHelper.h"

#include "KernelAlignment/KernelFactory.h"
#include "QuickConsistencyStage.h"

#include <thread>
#include <cmath>
#include <omp.h>


using namespace std;
using namespace quickprobs;

/// <summary>
/// See declaration for all the details.
/// </summary>
QuickConsistencyStage::QuickConsistencyStage(
	std::shared_ptr<clex::OpenCL> cl, 
	std::shared_ptr<Configuration> config) : ConsistencyStage(config), cl(cl)
{
	std::vector<std::string> defines;
	std::vector<int> files;

	files.push_back(KernelRepository::Utils);
	files.push_back(KernelRepository::MemoryTransfers);
	files.push_back(KernelRepository::ScoreType);
	files.push_back(KernelRepository::AuxiliaryStructures);
	files.push_back(KernelRepository::SparseMatrix);
	files.push_back(KernelRepository::Random);
	files.push_back(KernelRepository::MultipleConsistencyUtils);
	files.push_back(KernelRepository::MultipleConsistency_old);

	// adjust workgroup size if necessary
	if (config->stripeCount * config->stripeLength > cl->mainDevice->info->maxWorkgroupSize) {
		config->stripeCount = cl->mainDevice->info->maxWorkgroupSize / config->stripeLength;
	}

	defines.push_back("POSTERIOR_CUTOFF=" + std::to_string(config->algorithm.posteriorCutoff));
	defines.push_back("BINS_COUNT=" + std::to_string(config->binsCount));
	
	defines.push_back("STRIPE_COUNT=" + std::to_string(8));
	defines.push_back("STRIPE_LENGTH=" + std::to_string(8));
	defines.push_back("STRIPE_LENGTH_LOG2=" + std::to_string(3));

	Configuration::Algorithm::Consistency& cnf = config->algorithm.consistency;


	if (cnf.function == SelectivityFunction::Sum) {
		defines.push_back("SELECTIVITY_SUM");
	} else if (cnf.function == SelectivityFunction::Min) {
		defines.push_back("SELECTIVITY_MIN");
	} else if (cnf.function == SelectivityFunction::Max) {
		defines.push_back("SELECTIVITY_MAX");
	} else if (cnf.function == SelectivityFunction::Avg) {
		defines.push_back("SELECTIVITY_AVG");
	} 


	if (cnf.filter == SelectivityFilter::Deterministic) {
		defines.push_back("FILTER_VERTICAL");
	} else if (cnf.filter == SelectivityFilter::TriangleLowpass || cnf.filter == SelectivityFilter::TriangleHighpass) {
		defines.push_back("FILTER_LINEAR");
	} else if (cnf.filter == SelectivityFilter::HomographLowpass) {
		defines.push_back("FILTER_HOMOGRAPHIC");
	} else if (cnf.filter == SelectivityFilter::TriangleMidpass) {
		defines.push_back("FILTER_TRIANGLE_MIDPASS");
	}

	kernel_postprocessUnfiltered = KernelFactory::instance(cl).create(files, "MultipleConsistency_Postprocess", defines);
	kernel_postprocessFiltered = KernelFactory::instance(cl).create(files, "MultipleConsistency_PostprocessFiltered", defines);

	relaxKernels.push_back(std::shared_ptr<RelaxationKernelWrapper>(new RelaxationKernelWrapper(
		cl, "MultipleConsistency_RelaxOld", files, defines, 8, 8, 0)));

}

/// <summary>
/// See declaration for all the details.
/// </summary>
void QuickConsistencyStage::run(
	const float* seqsWeights,
	quickprobs::ISequenceSet& set,  
	Array<float>& distances, 
	Array<SparseMatrixType*>& sparseMatrices)
{
	STATS_WRITE("time.3-1-task preparation", 0.0);
	STATS_WRITE("time.3-2-relaxation", 0.0);
	STATS_WRITE("time.3-3-sparse unpacking", 0.0);

	//calculate sequence pairs for openmp model
	omp_set_num_threads(config->hardware.numThreads);
	
	auto& cms = dynamic_cast<ContiguousMultiSequence&>(set);
	const int numSeqs = cms.count();
	
	for (int iter = 0; iter < iterations; ++iter) {
		bool filterFlag = true; 
		
		// turn of filtering for last iteration if num filterings below 0
		if (config->algorithm.consistency.numFilterings < 0) {
			if (iter == iterations - 1) {
				filterFlag = false;
			}
		} else {
			filterFlag = iter < config->algorithm.consistency.numFilterings;
		}
		
		doRelaxationGpu(seqsWeights, cms, distances, sparseMatrices, filterFlag);
	}
}

/// <summary>
/// See declaration for all the details.
/// </summary>
void QuickConsistencyStage::doRelaxationGpu(
	const float* seqsWeights, 
	quickprobs::ContiguousMultiSequence& cms, 
	Array<float>& distances,
	Array<SparseMatrixType*>& sparseMatrices,
	bool filter)
{
	double timeUnpacking = 0;
	double timeKernel = 0;
	
//	printDistanceHistogram(distances);

	std::vector<float> weights(seqsWeights, seqsWeights + cms.count());

	TIMER_CREATE(prepareTimer);
	TIMER_CREATE(relaxTimer);
	TIMER_START(prepareTimer);
	
	int numSeqs = cms.count();

	// change weight vector to distance matrix
	int weightCount = numSeqs;
	
	std::vector<float> distArray = distances.getData();

	distArray[0 * numSeqs + 0] = selectivity.filter_a;
	distArray[1 * numSeqs + 1] = selectivity.filter_b;
	
	// selfweight in last element
	distArray[numSeqs * numSeqs - 1] = selfweight;

	cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
	clex::Buffer buffer_lengths(*cl, flags, numSeqs * sizeof(unsigned int), cms.lengths.data(), "buffer_lengths");
	clex::Buffer buffer_sortingMap(*cl, flags, numSeqs * sizeof(unsigned int), cms.sortingMap.data(), "buffer_sortingMap");
	clex::Buffer buffer_seqsWeights(*cl, flags, weightCount * sizeof(float), weights.data(), "buffer_seqsWeights");
	clex::Buffer buffer_distances(*cl, flags, distArray.size() * sizeof(float), distArray.data(), "buffer_distArray");
	clex::Buffer buffer_sparseOffsets(*cl, CL_MEM_READ_ONLY, numSeqs * numSeqs * sizeof(unsigned int), NULL, "buffer_sparseOffsets");

	// prepare waves
	std::vector<unsigned int> sparseOffsets;
	::size_t maxSectorBytes;
	auto sectors = generateSectors(sparseMatrices, cms, distances, sparseOffsets, maxSectorBytes);
	
	//printSelectivityHistogram(sectors);
	
	int sectorResolution = sqrt(sectors.size());

	auto maxTasksSector = *max_element(sectors.begin(), sectors.end(),
		[](std::shared_ptr<RelaxationSector>& w1, std::shared_ptr<RelaxationSector>& w2) -> bool {
			return w1->size() < w2->size();
	});

	LOG_DEBUG 
		<< endl << "CONSISTENCY MEMORY REPORT:" << endl
		<< "allocated MB: " << MemoryTools::processCurrentVirtual() / 1e6 << endl;

	clex::Buffer buffer_tasks(*cl, CL_MEM_READ_WRITE, maxTasksSector->size() * sizeof(RelaxationTask), NULL, "buffer_tasks");
	clex::Buffer buffer_sector_Sxy(*cl, CL_MEM_READ_WRITE, maxSectorBytes, NULL, "buffer_sector_Sxy");
	cl->mainDevice->mainQueue->enqueueWriteBuffer(buffer_sparseOffsets, CL_FALSE, 0, buffer_sparseOffsets.getBytesRequested(), sparseOffsets.data());

	std::shared_ptr<clex::Buffer> buffers_sector_Sxz[2];
	std::shared_ptr<clex::Buffer> buffers_sector_Szy[2];
	int streamsCount = std::min(sectorResolution, 2);

	for (int i = 0; i < streamsCount; ++i) {
		buffers_sector_Sxz[i] = std::shared_ptr<clex::Buffer>(new clex::Buffer(*cl, CL_MEM_READ_ONLY, maxSectorBytes, NULL, "buffer_sector_Sxz"));
		buffers_sector_Szy[i] = (sectorResolution > 1) ?
			std::shared_ptr<clex::Buffer>(new clex::Buffer(*cl, CL_MEM_READ_ONLY, maxSectorBytes, NULL, "buffer_sector_Szy")) :
			buffers_sector_Sxz[i];
	}

	LOG_DEBUG 
		<< "dynamic CL buffers MB: " << clex::Buffer::getTotalBytesAllocated() / 1e6 << endl
		<< "allocated MB: " << MemoryTools::processCurrentVirtual() / 1e6 << endl;

	// set parameters for all relax kernels
	for (auto& kernel : relaxKernels) {
		clCall(kernel->obj->setArg(0, buffer_lengths));
		clCall(kernel->obj->setArg(1, buffer_sortingMap));
		clCall(kernel->obj->setArg(2, buffer_seqsWeights));
		clCall(kernel->obj->setArg(3, buffer_distances));
		clCall(kernel->obj->setArg(4, buffer_sparseOffsets));
		clCall(kernel->obj->setArg(5, buffer_tasks));
		clCall(kernel->obj->setArg(6, buffer_sector_Sxy));
	}

	std::shared_ptr<clex::Kernel> kernel_postprocess = filter ? kernel_postprocessFiltered : kernel_postprocessUnfiltered;

	clCall(kernel_postprocess->setArg(0, buffer_lengths));
	clCall(kernel_postprocess->setArg(1, buffer_sparseOffsets));
	clCall(kernel_postprocess->setArg(2, buffer_tasks));
	clCall(kernel_postprocess->setArg(3, buffer_sector_Sxy));
	clCall(kernel_postprocess->setArg(4, config->algorithm.posteriorCutoff));

	TIMER_STOP(prepareTimer);

	TIMER_START(relaxTimer);
	
	int sectorSize = mathex::ceildiv(numSeqs, sectorResolution);
	int kernelExecution = 0;
	int streamId = 0;
	
	std::vector<cl::Event> relaxFinished(streamsCount * config->binsCount);
	std::vector<cl::Event> copyFinished(streamsCount);
	cl::Event filterFinished;
	cl::Event preprocessFinished;
	std::thread unpackThread;

	buffer_type *outputRawData = new buffer_type[maxSectorBytes / sizeof(buffer_type)];

	std::vector<buffer_type*> pool;

	if (config->algorithm.consistency.copy == SectorCopy::Row) {
		pool.resize(sectorResolution + streamsCount);
	} else if (config->algorithm.consistency.copy == SectorCopy::None) {
		pool.resize(1 + 2 * streamsCount);
	}

	for (auto& p : pool) { p = new buffer_type[maxSectorBytes / sizeof(buffer_type)]; }

	LOG_DEBUG 
		<< "input host buffers MB: " << pool.size() * maxSectorBytes / 1e6 << endl
		<< "output host buffer MB: " << maxSectorBytes / 1e6 << endl
		<< "allocated MB: " << MemoryTools::processCurrentVirtual() / 1e6 << endl;

	for (int sectorId = 0; sectorId < sectors.size(); ++sectorId) {

		if (unpackThread.joinable()) {
			unpackThread.join();
		}
	
		int sector_i = sectorId / sectorResolution;
		int sector_j = sectorId % sectorResolution;
		auto& sector_Sxy = *sectors[sectorId];
		
		if (config->algorithm.consistency.copy == SectorCopy::None) {
			for (int k = 0; k < sectorResolution; ++k) { // reset horizontal and vertical
				sectors[sector_i * sectorResolution + k]->resetRawData();
				sectors[k * sectorResolution + sector_j]->resetRawData();
			}
			sector_Sxy.packExternal(sparseMatrices, sparseOffsets, cms.sortingMap, pool[0]); // does nothing if sector already packed
		
		} else if (config->algorithm.consistency.copy == SectorCopy::Row) {
			if (sector_j == 0) { // reset horizontal
				for (int k = 0; k < sectorResolution; ++k) {
					sectors[sectorId + k]->resetRawData();
				}
			}
			for (int k = 0; k < sectorResolution; ++k) { // reset vertical
				if (k != sector_i) {
					sectors[k * sectorResolution + sector_j]->resetRawData();
				}
			}
			sector_Sxy.packExternal(sparseMatrices, sparseOffsets, cms.sortingMap, pool[sector_j]); // does nothing if sector already packed
		
		} else {
			sector_Sxy.packInternal(sparseMatrices, sparseOffsets, cms.sortingMap);
		}

		// omit sector if no tasks defined
		if (sector_Sxy.getTasks().size() == 0) {
			continue;
		}

		// copy sector-specific data
		clCall(cl->mainDevice->mainQueue->enqueueWriteBuffer(buffer_sector_Sxy, CL_FALSE, 0, sector_Sxy.getElementsCount() * sizeof(buffer_type), sector_Sxy.getRawData())); // start transfer
		clCall(cl->mainDevice->mainQueue->enqueueWriteBuffer(buffer_tasks, CL_FALSE, 0, sector_Sxy.getTasks().size() * sizeof(RelaxationTask),  
			sector_Sxy.getTasks().data())); // start transfer

		LOG_DEBUG << endl << "RELAXATION SECTOR " << sectorId + 1 << "/" << sectors.size() << ": "
			<< "sparse_elems= " << sector_Sxy.getElementsCount()
			<< ", group_size= " << config->stripeLength * config->stripeCount
			<< ", max_sparse= " << sector_Sxy.maxSparseRow()
			<< ", max_dense_hor= " << sector_Sxy.maxHorizontal()
			<< ", max_dense_ver= " << sector_Sxy.maxVertical()
			<< ", stripe_count= " << config->stripeCount << endl;

		int start_k = 0;//max(0, sector_j - 1);

		for (int k = start_k; k < start_k + sectorResolution; ++k, ++kernelExecution) {
			
			int sequenceSubset = k % sectorResolution;
			
			auto &sector_Sxz = *sectors[sector_i * sectorResolution + sequenceSubset];
			auto &sector_Szy = *sectors[sequenceSubset * sectorResolution + sector_j];

			streamId = kernelExecution % streamsCount;

			// wait for current to finish
			if (kernelExecution > 1) {
				for (int r = 0; r < relaxKernels.size(); ++r) {
					relaxFinished[streamId * config->binsCount + r].wait();
				}
			}
			
			// copy data to first buffer
			if (config->algorithm.consistency.copy == SectorCopy::None) {
				sector_Sxz.packExternal(sparseMatrices, sparseOffsets, cms.sortingMap, pool[1 + streamId]);
			} else if (config->algorithm.consistency.copy == SectorCopy::Row) {
				sector_Sxz.packExternal(sparseMatrices, sparseOffsets, cms.sortingMap, pool[sequenceSubset]);
			} else {
				sector_Sxz.packInternal(sparseMatrices, sparseOffsets, cms.sortingMap);
			}

			clCall(cl->mainDevice->queues[streamId * config->binsCount]->enqueueWriteBuffer(*buffers_sector_Sxz[streamId], CL_FALSE, 0, 
				sector_Sxz.getElementsCount() * sizeof(buffer_type),  sector_Sxz.getRawData(), NULL, &copyFinished[0])); // start transfer

			// copy data to second buffer only when necessary
			if (streamsCount > 1) {
				
				if (config->algorithm.consistency.copy == SectorCopy::None) {
					sector_Szy.packExternal(sparseMatrices, sparseOffsets, cms.sortingMap, pool[3 + streamId]);
				} else if (config->algorithm.consistency.copy == SectorCopy::Row) {
					sector_Szy.packExternal(sparseMatrices, sparseOffsets, cms.sortingMap, pool[sectorResolution + streamId]);
				} else {
					sector_Szy.packInternal(sparseMatrices, sparseOffsets, cms.sortingMap);
				}
				clCall(cl->mainDevice->queues[streamId * config->binsCount]->enqueueWriteBuffer(*buffers_sector_Szy[streamId], CL_FALSE, 0, 
					sector_Szy.getElementsCount() * sizeof(buffer_type),  sector_Szy.getRawData(), NULL, &copyFinished[1])); // start transfer
			} 

			// set proper kernel arguments
			for (auto& kernel : relaxKernels) {
				clCall(kernel->obj->setArg(7, *buffers_sector_Sxz[streamId]));
				clCall(kernel->obj->setArg(8, *buffers_sector_Szy[streamId]));
			}

			// relaxation kernel
			RelaxationSchedule relaxSchedule;
			relaxSchedule.numSeqs = numSeqs;
			relaxSchedule.sparseWidth = sector_Sxy.maxSparseRow() + 1;
			relaxSchedule.sequenceBegin = sequenceSubset * sectorSize;
			relaxSchedule.sequenceEnd = std::min(relaxSchedule.sequenceBegin + sectorSize, numSeqs);

			cl::NDRange globalRange(config->stripeLength * config->stripeCount * sector_Sxy.getTasks().size());
			cl::NDRange localRange(config->stripeLength * config->stripeCount);

			int localBytes = config->stripeCount * config->stripeLength * sizeof(OpenCLSparseMatrix::cell_type) +
				config->stripeCount * relaxSchedule.sparseWidth * sizeof(OpenCLSparseMatrix::cell_type);
			clCall(relaxKernels[0]->obj->setArg(9, relaxSchedule));
			clCall(relaxKernels[0]->obj->setArg(10, cl::__local(localBytes)));

			LOG_DEBUG << "Relax old " << sequenceSubset + 1 << ": "  << "local_mem= " << relaxKernels[0]->obj->localMemSize << "...";
			clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
				*relaxKernels[0]->obj, cl::NDRange(0), globalRange, localRange, &copyFinished, &relaxFinished[streamId * config->binsCount]));
			LOG_DEBUG << "scheduled..." << endl;
		}
		
		for (int s= 0; s < streamsCount; ++s) {
			for (int r = 0; r < relaxKernels.size(); ++r) {
				relaxFinished[s * config->binsCount + r].wait();
			}
		}

		RelaxationSchedule schedule;
		schedule.numSeqs = numSeqs;
		schedule.sparseWidth = sector_Sxy.maxSparseRow() + 1;
		schedule.sequenceBegin = sector_Sxy.maxLength() + 1;
	
		int localFilteringBytes = filter 
			? config->stripeCount * schedule.sparseWidth * sizeof(OpenCLSparseMatrix::cell_type)
			: 4;

		// wait for previous unpack to finish if necessary
		if (unpackThread.joinable()) {
			unpackThread.join();
		}
		
		LOG_DEBUG << "Filtering: " << "local_mem= " << localFilteringBytes << "...";
		clCall(kernel_postprocess->setArg(5, schedule));
		clCall(kernel_postprocess->setArg(6, cl::__local(localFilteringBytes)));

		clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
			*kernel_postprocess, 
			cl::NDRange(0), 
			cl::NDRange(sector_Sxy.getTasks().size() * config->stripeLength * config->stripeCount), 
			cl::NDRange(config->stripeLength * config->stripeCount)));
		
		clCall(cl->mainDevice->mainQueue->enqueueReadBuffer(
			buffer_sector_Sxy,
			CL_FALSE,
			0,
			sector_Sxy.getElementsCount() * sizeof(buffer_type),
			outputRawData,
			NULL,
			&filterFinished));
		LOG_DEBUG << "scheduled..." << endl;

		// run sector unpack thread
		unpackThread = std::thread([&filterFinished, &sector_Sxy, &sparseOffsets, outputRawData, &sparseMatrices, &timeUnpacking](){
			LOG_DEBUG << "Unpacking thread started. Waiting for filtering kernel to finish...";
			filterFinished.wait();
			LOG_DEBUG << "filtered" << endl << "Unpacking sector...";
			TIMER_CREATE(timer);
			TIMER_START(timer);
			sector_Sxy.unpack(outputRawData, sparseOffsets, sparseMatrices);
			TIMER_STOP(timer);
			timeUnpacking += timer.seconds();
			LOG_DEBUG << "unpacked" << endl;

		}); 
	}

	TIMER_STOP(relaxTimer);

	unpackThread.join();

	STATS_WRITE("scheduling.relaxation waves", sectors.size());
	STATS_ADD("time.3-1-task preparation", prepareTimer.seconds());
	STATS_ADD("time.3-2-relaxation", relaxTimer.seconds());
	STATS_ADD("time.3-3-sparse unpacking", timeUnpacking);
	
	for (auto& p : pool) { delete [] p; }
	delete [] outputRawData;
}



std::vector<std::shared_ptr<RelaxationSector>> quickprobs::QuickConsistencyStage::generateSectors(
	const Array<SparseMatrixType*>& sparseMatrices, 
	const ContiguousMultiSequence& sequences,
	const Array<float>& distances,
	std::vector<unsigned int>& sparseOffsets,
	::size_t& maxSectorBytes)
{
	std::vector<std::shared_ptr<RelaxationSector>> sectors;

	// calculate sparse memory requirements
	::size_t openClSparseMemory = 0;
	int i = 0;
	#pragma omp parallel for private(i) default(shared) reduction(+ : openClSparseMemory) schedule(static) 
	for (i = 0; i < sequences.count(); ++i) {
		for (int j = 0; j < sequences.count(); ++j) {	
			if (i != j) {
				auto Sij = sparseMatrices[i][j];
				openClSparseMemory += OpenCLSparseMatrix::bytesNeeded(Sij->getSeq1Length(), Sij->getNumCells());
			}
		}
	}

	// sector must fit in sixth of total memory size and half of max alloc size
	::size_t metaBytes = cl->mainDevice->info->alignUp(sequences.count() * sequences.count() * sizeof(RelaxationTask));
	::size_t maxAllocBytes = config->hardware.gpuMemFactor * cl->mainDevice->info->maxAllocSize_Corrected;
	::size_t availableBytes = config->hardware.gpuMemFactor * cl->mainDevice->info->globalMemSize - (clex::Buffer::getTotalBytesAllocated() + metaBytes);
	::size_t allowedBytes = std::min(availableBytes / 10, maxAllocBytes / 2); 
	openClSparseMemory = mathex::ceilround(openClSparseMemory, (size_t)cl->mainDevice->info->memAddrAlign);
	int sectorsCount = mathex::ceildiv(openClSparseMemory, allowedBytes);

	LOG_DEBUG << endl << "GENERATING SECTORS: " << endl
		<< "opencl sparse MB = " << openClSparseMemory / 1e6 << endl
		<< "initial allowed sector MB = " << allowedBytes / 1e6 << endl;
		
	// update allowed bytes
	allowedBytes = cl->mainDevice->info->alignDown(availableBytes / 6); 

	LOG_DEBUG << "final allowed sector MB = " << allowedBytes / 1e6 << endl;

	// fixme: ugly hardcoded tuning
	int resolution = 0;
	
//	if (sectorsCount < 10) {
//		sectorResolution = sectorsCount;
//	} else if (sectorsCount < 100) {
//		sectorResolution = 10;
//	} else {
		
//	}

	std::shared_ptr<RelaxationSector> maxMemorySector; 
	int iteration = 1;

	do {
		resolution = (::size_t)ceil(sqrt((float)sectorsCount));
		sectors = RelaxationSector::generate(sparseMatrices, sequences.sortingMap, distances, selectivity, sparseOffsets, resolution);
		maxMemorySector = *max_element(sectors.begin(), sectors.end(),
			[](std::shared_ptr<RelaxationSector>& w1, std::shared_ptr<RelaxationSector>& w2) -> bool {
				return w1->getElementsCount() < w2->getElementsCount();
		});
		
		maxSectorBytes = mathex::ceilround(maxMemorySector->getBytes(), (size_t)cl->mainDevice->info->memAddrAlign);
		LOG_DEBUG 
			<< "iteration " << iteration << ":" << endl
			<< "current max sector MB = " << maxSectorBytes / 1e6 << endl
			<< "resolution = " << resolution << endl;

		sectorsCount = sectorsCount * 2;
		++iteration;

	} while (maxSectorBytes > allowedBytes);

	return sectors;
}

void quickprobs::QuickConsistencyStage::printSelectivityHistogram(
	std::vector<std::shared_ptr<RelaxationSector>> sectors)
{
	std::vector<std::pair<int,int>> selectivityHisto;
	selectivityHisto.push_back(std::pair<int,int>(0,0));
	selectivityHisto.push_back(std::pair<int,int>(1,0));
	selectivityHisto.push_back(std::pair<int,int>(2,0));
	selectivityHisto.push_back(std::pair<int,int>(3,0));
	selectivityHisto.push_back(std::pair<int,int>(5,0));
	selectivityHisto.push_back(std::pair<int,int>(10,0));

	
	for (int i = 20; i <= config->algorithm.consistency.selectivity; i += 20) {
		selectivityHisto.push_back(std::pair<int,int>(i, 0));
	}
	selectivityHisto.push_back(make_pair<int, int>(100000, 0));
	int attempts = 0;

	for (const auto& s: sectors) {
		for (const auto& t: s->getTasks()) {	
			auto it = find_if(selectivityHisto.begin(), selectivityHisto.end(), [&t](std::pair<int,int>& bin)->bool {
				return t.acceptedCount <= bin.first;
			});	

			++it->second;
			++attempts;
		}
	}

	LOG_DEBUG << endl << "Selectivity histogram(" <<  attempts << ")" << endl;
	for (auto& bin: selectivityHisto) {
		LOG_DEBUG << bin.first << ": " << bin.second << endl;
	}

	getchar();
}

void quickprobs::QuickConsistencyStage::printDistanceHistogram(Array<float>& distances)
{
	int numSeqs = distances.size();
	std::vector<std::pair<float,int>> distanceHisto(20);
	std::vector<std::pair<float,int>> pairsHisto(20);
	//float maxVal = std::numeric_limits<float>::max();
	float maxVal = (float)numSeqs;

	float current = 10.0f;
	float step = 10.0f;
	for (auto &x : distanceHisto) {
		x.first = current;
		current += step;
	}
	distanceHisto.push_back(make_pair<float, int>(100000, 0));

	current = 10.0f;
	for (auto &x : pairsHisto) {
		x.first = current;
		current += step;
	}
	pairsHisto.push_back(make_pair<float,int>(100000, 0));

	int distCount = 0;
	int pairsCount = 0;

	for (int i = 0; i < distances.size(); ++i) {
		for (int j = i + 1; j < distances.size(); ++j) {
		
			for (int k = 0; k < distances.size(); ++k) {
				float sum = distances[i][k] + distances[k][j];
				auto it = find_if(pairsHisto.begin(), pairsHisto.end(), [sum](std::pair<float,int>& bin)->bool {
					return sum < bin.first;
				});
				++it->second;
				++pairsCount;
			}

			float d = distances[i][j]; 
			auto it = find_if(distanceHisto.begin(), distanceHisto.end(), [d](std::pair<float,int>& bin)->bool {
				return d < bin.first;
			});
			++it->second;
			++distCount;
		}
	}

	LOG_DEBUG << endl << "Sequences count: " << numSeqs << endl;

	LOG_DEBUG << endl << "Distance histogram (" <<  distCount << ")" << endl;
	for (auto& bin: distanceHisto) {
		LOG_DEBUG << bin.first << ": " << bin.second << endl;
	}

	LOG_DEBUG << endl << "Distance sums histogram (" <<  pairsCount << ")" << endl;
	for (auto& bin: pairsHisto) {
		LOG_DEBUG << bin.first << ": " << bin.second << endl;
	}

	getchar();
}

void quickprobs::QuickConsistencyStage::printSparseRowsHistogram(Array<SparseMatrixType*>& sparseMatrices)
{
	std::map<int,int> histogram;
	histogram[0] = 0;
	histogram[1] = 0;
	histogram[2] = 0;
	histogram[4] = 0;
	histogram[8] = 0;
	histogram[16] = 0;
	histogram[32] = 0;
	SparseHelper::getHistogram(sparseMatrices, histogram);
	cout << "histogram:" << endl;
	for (auto& h : histogram) {
		cout << "h[" << h.first << "] = " << h.second << endl; 
	}
}
