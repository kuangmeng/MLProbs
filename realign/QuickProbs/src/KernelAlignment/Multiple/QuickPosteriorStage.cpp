#define _VARIADIC_MAX 6

#include "Common/mathex.h"
#include "Common/Log.h"
#include "Common/MemoryTools.h"

#include "Alignment/DataStructures/ISequenceSet.h"
#include "Alignment/DataStructures/Sequence.h"
#include "Alignment/DataStructures/MultiSequence.h"
#include "Alignment/DataStructures/ContiguousMultiSequence.h"
#include "Alignment/DataStructures/SparseMatrixType.h"
#include "Alignment/DataStructures/SparseHelper.h"
#include "Alignment/Pairwise/ProteinHmm5.h"
#include "Alignment/Multiple/ProbabilisticParams.h"
#include "Alignment/Multiple/ParallelProbabilisticModel.h"
#include "Alignment/Multiple/PartitionFunction.h"
#include "Alignment/Multiple/Configuration.h"
#include "Alignment/Multiple/ScoreType.h"

#include "QuickPosteriorStage.h"
#include "PosteriorTasksWave.h"

#include <omp.h>

#include <thread>
#include <atomic>

using namespace std;
using namespace quickprobs;

QuickPosteriorStage::QuickPosteriorStage(
	std::shared_ptr<clex::OpenCL> cl,
	std::shared_ptr<Configuration> config) : PosteriorStage(config), cl(cl), kernelSets(2), maxWorkgroupSizes(2)
{

	kernelSets[0].fillAll(*config, cl);
	kernelSets[1].fillAllLong(*config, cl);

	size_t localKernelMem = 0;
	for (auto& set: kernelSets) {
		for (auto& k: set.getItems())
			if (k != nullptr && k->localMemSize > localKernelMem) {
				localKernelMem = k->localMemSize;
		}
	}

	// get max workgroup size for classic kernels
	size_t localMemAvailable = cl->mainDevice->info->localMemSize - localKernelMem;
	size_t bytesPerSymbol = PosteriorTask::elementsPerSymbol(config->optimisation.useDoublePartition) * sizeof(buffer_type);

	this->maxWorkgroupSizes[0] = std::min(kernelSets[0].maxWorkgroupSize(), localMemAvailable / bytesPerSymbol);
	this->maxWorkgroupSizes[1] = std::min(kernelSets[1].maxWorkgroupSize(), localMemAvailable / bytesPerSymbol);
}

void QuickPosteriorStage::run(
	quickprobs::ISequenceSet& set,
	Array<float>& distances,
	Array<SparseMatrixType*>& matrices)
{
	omp_set_num_threads(config->hardware.numThreads);
	LOG_NORMAL << "Posterior with OpenMP using " << omp_get_max_threads() << " threads..." << endl;

	FilteredSparseMatrix::initialise(config->algorithm.posteriorCutoff, 1024);
	auto& cms = dynamic_cast<ContiguousMultiSequence&>(set);
	const int numSeqs = cms.count();

	int maxSparseRow = RelaxationTask::countMaxSparseLength(
		cl->mainDevice->info->localMemSize - 128,  //FIXME - assumed local memory requirement for relaxation kernel
		config->stripeCount,
		config->stripeLength);
	PackedSparseMatrix::setSparseRowThreshold(maxSparseRow);


	::size_t wavesCount[3] = {0};
	::size_t tasksCount[3] = {0};
	std::atomic< ::size_t> sparseBytes[3];
	sparseBytes[0] = sparseBytes[1] = sparseBytes[2] = 0;

	// prepare data for passing them to kernel
	timeCpuProcessing = 0;
	TIMER_CREATE(prepareTimer);
	TIMER_CREATE(denseTimer);
	TIMER_CREATE(sparseTimer);
	TIMER_CREATE(unpackTimer);
	TIMER_CREATE(auxTimer);

	TIMER_START(prepareTimer);

	LOG_DEBUG << endl << "POSTERIOR MEMORY REPORT:" << endl
		<< "allocated MB: " << MemoryTools::processCurrentVirtual() / 1e6 << endl;

	// allocate wave-common buffers
	cl_mem_flags flags = CL_MEM_READ_ONLY | CL_MEM_USE_HOST_PTR;
	clex::Buffer buffer_sequences(*cl, flags, cms.sequenceData.size() * sizeof(char), cms.sequenceData.data(), "buffer_sequences");
	clex::Buffer buffer_offsets(*cl, flags, cms.offsets.size() * sizeof(unsigned int), cms.offsets.data(), "buffer_offsets");
	clex::Buffer buffer_lengths(*cl, flags, cms.lengths.size() * sizeof(unsigned int), cms.lengths.data(), "buffer_lengths");
	clex::Buffer buffer_probParams(*cl, flags, model->params.sizeInBytes(), model->params.data(), "buffer_probParams");
	clex::Buffer buffer_funcParams(*cl, flags, function->params->sizeInBytes(), function->params->data(), "buffer_funcParams");

	LOG_DEBUG
		<< "static CL buffers MB: " << clex::Buffer::getTotalBytesAllocated() / 1e6 << endl;

	// generate waves
	::size_t availableGlobalMemory = (::size_t)(cl->mainDevice->info->globalMemSize - clex::Buffer::getTotalBytesAllocated());
	::size_t maxBufferSize = cl->mainDevice->info->maxAllocSize_Corrected - 16 * cl->mainDevice->info->memAddrAlign;

	auto tasks = std::move(PosteriorTasksWave::generateTasks(
		cms,
		distances,
		availableGlobalMemory,
		maxBufferSize,
		this->maxWorkgroupSizes[0]));

	int columnElements = (tasks[1].size() > 0) ? ((cms.maxLength + 1) * 6) : 0;

	auto waves = std::move(PosteriorTasksWave::generate(
		tasks[0],
		availableGlobalMemory * 0.9,
		std::min(std::min((size_t)512e6, maxBufferSize), availableGlobalMemory / 4),
		columnElements,
		0));

	auto longWaves = std::move(PosteriorTasksWave::generate(
		tasks[1],
		availableGlobalMemory * 0.9,
		std::min(maxBufferSize, availableGlobalMemory / 4),
		columnElements,
		1));

	auto veryLongWaves = std::move(PosteriorTasksWave::generate(
		tasks[2],
		std::numeric_limits< ::size_t>::max(),
		std::numeric_limits< ::size_t>::max(),
		columnElements,
		2));

	wavesCount[0] = waves.size();
	wavesCount[1] = longWaves.size();
	wavesCount[2] = veryLongWaves.size();

	std::thread cpuWaveThread;
	if (veryLongWaves.size() > 0) {
		auto wave = *(veryLongWaves.begin());
		tasksCount[2] = wave->size();

		cpuWaveThread = std::thread(
			&QuickPosteriorStage::computeWaveCpu,
			this,
			std::ref(*wave),
			std::ref(cms),
			std::ref(distances),
			std::ref(matrices),
			std::ref(sparseBytes[2]));
	}

	waves.insert(waves.end(), longWaves.begin(), longWaves.end());

	// get some wave statistics
	auto waveMaxMemory = waves.front();
	auto waveMaxTasks = waves.front();
	int64_t minSparseBytes = 0;
	denseBytes = 0;
	for (auto& w: waves ) {
		if (waveMaxMemory->posteriorElemsCount < w->posteriorElemsCount) { waveMaxMemory = w; }
		if (waveMaxTasks->size() < w->size()) { waveMaxTasks = w; }
		minSparseBytes += w->minSparseBytes;
		denseBytes += 2 * w->posteriorElemsCount * sizeof(buffer_type);
	}
	// fixme: sparse matrix allocation overhead
	double overhead = 1.15;
	int64_t reservedBytes = static_cast<int64_t>(overhead * minSparseBytes + cl->mainDevice->info->globalMemSize);
	int64_t extraSparseBytes = static_cast<int64_t>(config->hardware.memoryLimitMb * 1e6) - reservedBytes;

	double fillRate = (double)(extraSparseBytes) / (this->denseBytes - minSparseBytes);
	fillRate = std::min(std::max(0.0, fillRate / overhead), 1.0);

	// allocate buffers depending on a largest wave
	int kernelsCount = 1;

	// round size to the multiplicity of 2 * memAddrAlign
	size_t auxiliaryElementsCount =  mathex::ceilround(waveMaxMemory->auxiliaryElemsCount,
		(size_t)cl->mainDevice->info->memAddrAlign * 2 / sizeof(buffer_type));

	std::shared_ptr<clex::Buffer> buffer_columns;
	clex::Buffer buffer_auxiliary(*cl, CL_MEM_READ_WRITE, auxiliaryElementsCount * sizeof(buffer_type), NULL, "buffer_auxiliaryMatrices");
	clex::Buffer buffer_posterior(*cl, CL_MEM_READ_WRITE, waveMaxMemory->posteriorElemsCount * sizeof(buffer_type), NULL, "buffer_posteriorMatrices");
	clex::Buffer buffer_tasks(*cl, CL_MEM_READ_WRITE, waveMaxTasks->size() * sizeof(PosteriorTask), NULL, "buffer_tasks");
	clex::Buffer buffer_tasksSparse(*cl, CL_MEM_READ_WRITE, waveMaxTasks->size() * sizeof(PosteriorTask), NULL, "buffer_tasksSparse");

	LOG_DEBUG
		<< "dense MB = " << denseBytes / 1e6 << endl
		<< "memory limit MB = " << config->hardware.memoryLimitMb << endl
		<< "minimum sparse MB = " << minSparseBytes / 1e6 << endl
		<< "reserved (sparse + overhead + consistency buffers) MB = " << reservedBytes / 1e6 << endl
		<< "extra sparse MB = " << extraSparseBytes / 1e6 << endl
		<< "fill rate (extra sparse - overhead) = " << fillRate << endl
		<< "dynamic CL buffers MB: " << clex::Buffer::getTotalBytesAllocated() / 1e6 << endl
		<< "auxiliary CL buffer MB = " << buffer_auxiliary.getBytesAllocated() / 1e6 << endl
		<< "posterior CL buffer MB = " << buffer_posterior.getBytesAllocated() / 1e6 << endl
		<< "allocated MB = " << MemoryTools::processCurrentVirtual() / 1e6 << endl
		<< "max short workgroup = " << maxWorkgroupSizes[0] << endl
		<< "max long workgroup = " << maxWorkgroupSizes[1] << endl << endl;

	// create subbuffer
	cl_buffer_region region = {buffer_auxiliary.getBytesAllocated() / 2, buffer_auxiliary.getBytesAllocated() / 2};
	cl::Buffer subbuffer_auxiliaryHi = buffer_auxiliary.createSubBuffer(CL_MEM_READ_WRITE, CL_BUFFER_CREATE_TYPE_REGION, &region);

	// if there are long kernels
	if (longWaves.size() > 0) {
		kernelsCount = 2;
		buffer_columns = std::shared_ptr<clex::Buffer>(new clex::Buffer(
			*cl, CL_MEM_READ_WRITE, waveMaxTasks->size() * columnElements * sizeof(buffer_type), NULL, "buffer_columns"));
	}

	for (int i = 0; i < kernelsCount; ++i) {
		auto& kernels = kernelSets[i];
		kernels(KernelSet::PARTITION_FORWARD).setArg(0, buffer_sequences);
		kernels(KernelSet::PARTITION_FORWARD).setArg(1, buffer_offsets);
		kernels(KernelSet::PARTITION_FORWARD).setArg(2, buffer_lengths);
		kernels(KernelSet::PARTITION_FORWARD).setArg(3, buffer_funcParams);
		kernels(KernelSet::PARTITION_FORWARD).setArg(4, buffer_auxiliary);
		kernels(KernelSet::PARTITION_FORWARD).setArg(5, buffer_posterior);
		kernels(KernelSet::PARTITION_FORWARD).setArg(6, buffer_tasks);

		kernels(KernelSet::PARTITION_REVERSE).setArg(0, buffer_sequences);
		kernels(KernelSet::PARTITION_REVERSE).setArg(1, buffer_offsets);
		kernels(KernelSet::PARTITION_REVERSE).setArg(2, buffer_lengths);
		kernels(KernelSet::PARTITION_REVERSE).setArg(3, buffer_funcParams);
		kernels(KernelSet::PARTITION_REVERSE).setArg(4, buffer_auxiliary);
		kernels(KernelSet::PARTITION_REVERSE).setArg(5, buffer_posterior);
		kernels(KernelSet::PARTITION_REVERSE).setArg(6, buffer_tasks);

		if (config->optimisation.divideHmmKernels) {
			kernels(KernelSet::HMM_FORWARD).setArg(0, buffer_sequences);
			kernels(KernelSet::HMM_FORWARD).setArg(1, buffer_offsets);
			kernels(KernelSet::HMM_FORWARD).setArg(2, buffer_lengths);
			kernels(KernelSet::HMM_FORWARD).setArg(3, buffer_probParams);
			kernels(KernelSet::HMM_FORWARD).setArg(4, buffer_auxiliary);
			kernels(KernelSet::HMM_FORWARD).setArg(5, buffer_posterior);
			kernels(KernelSet::HMM_FORWARD).setArg(6, buffer_tasks);

			kernels(KernelSet::HMM_BACKWARD).setArg(0, buffer_sequences);
			kernels(KernelSet::HMM_BACKWARD).setArg(1, buffer_offsets);
			kernels(KernelSet::HMM_BACKWARD).setArg(2, buffer_lengths);
			kernels(KernelSet::HMM_BACKWARD).setArg(3, buffer_probParams);
			kernels(KernelSet::HMM_BACKWARD).setArg(4, buffer_auxiliary);
			kernels(KernelSet::HMM_BACKWARD).setArg(5, buffer_posterior);
			kernels(KernelSet::HMM_BACKWARD).setArg(6, buffer_tasks);

			kernels(KernelSet::HMM_COMBINE).setArg(0, buffer_lengths);
			kernels(KernelSet::HMM_COMBINE).setArg(1, buffer_probParams);
			kernels(KernelSet::HMM_COMBINE).setArg(2, buffer_auxiliary);
			kernels(KernelSet::HMM_COMBINE).setArg(3, buffer_tasks);
		} else {
			kernels(KernelSet::HMM_ALL).setArg(0, buffer_sequences);
			kernels(KernelSet::HMM_ALL).setArg(1, buffer_offsets);
			kernels(KernelSet::HMM_ALL).setArg(2, buffer_lengths);
			kernels(KernelSet::HMM_ALL).setArg(3, buffer_probParams);
			kernels(KernelSet::HMM_ALL).setArg(4, buffer_auxiliary);
			kernels(KernelSet::HMM_ALL).setArg(5, buffer_posterior);
			kernels(KernelSet::HMM_ALL).setArg(6, buffer_tasks);
		}

		kernels(KernelSet::FINALIZATION).setArg(0, buffer_lengths);
		kernels(KernelSet::FINALIZATION).setArg(1, buffer_auxiliary);
		kernels(KernelSet::FINALIZATION).setArg(2, buffer_posterior);
		kernels(KernelSet::FINALIZATION).setArg(3, buffer_tasks);

		kernels(KernelSet::SPARSE).setArg(0, buffer_lengths);
		kernels(KernelSet::SPARSE).setArg(2, buffer_posterior);
		kernels(KernelSet::SPARSE).setArg(3, buffer_tasksSparse);

		// set column buffer if necessary
		if (i == 1) {
			kernels(KernelSet::PARTITION_FORWARD).setArg(9, *buffer_columns);
			kernels(KernelSet::PARTITION_REVERSE).setArg(9, *buffer_columns);
			if (config->optimisation.divideHmmKernels) {
				kernels(KernelSet::HMM_FORWARD).setArg(9, *buffer_columns);
				kernels(KernelSet::HMM_BACKWARD).setArg(9, *buffer_columns);
			} else {
				kernels(KernelSet::HMM_ALL).setArg(9, *buffer_columns);
			}
			kernels(KernelSet::FINALIZATION).setArg(6, *buffer_columns);

			kernels(KernelSet::PARTITION_FORWARD).setArg(10, columnElements);
			kernels(KernelSet::PARTITION_REVERSE).setArg(10, columnElements);
			if (config->optimisation.divideHmmKernels) {
				kernels(KernelSet::HMM_FORWARD).setArg(10, columnElements);
				kernels(KernelSet::HMM_BACKWARD).setArg(10, columnElements);
			} else {
				kernels(KernelSet::HMM_ALL).setArg(10, columnElements);
			}
			kernels(KernelSet::FINALIZATION).setArg(7, columnElements);
		}
	}

	buffer_type* auxiliaryMatrices = new buffer_type[auxiliaryElementsCount];

	LOG_DEBUG
		<< "auxiliary host buffer MB: " << auxiliaryElementsCount * sizeof(buffer_type) / 1e6 << endl
		<< "allocated MB: " << MemoryTools::processCurrentVirtual() / 1e6 << endl << endl;

	TIMER_STOP(prepareTimer);

	// calculate other waves on the GPU
	std::vector<cl::Event> finishedDense(1);
	std::vector<cl::Event> finishedSparse(2);
	std::vector<std::thread> unpackThreads(2);

	LOG_DEBUG << endl << "POSTERIOR STATISTICS:" << endl
		<< "short = " << wavesCount[0] << " waves, " << tasksCount[0] << " tasks" << endl
		<< "long = " << wavesCount[1] << " waves, " << tasksCount[1] << " tasks" << endl
		<< "very long = " << wavesCount[2] << " waves, " << tasksCount[2] << " tasks" << endl << endl;



	int waveId = 1;
	double timePhases[6] = {0};

	std::vector<cl::Event*> phaseFinished(6, nullptr);
	if (config->optimisation.kernelProfiling) {
		std::for_each(phaseFinished.begin(), phaseFinished.end(), [](cl::Event*& e) { e = new cl::Event(); });
	}

	for (auto ptr = waves.begin(); ptr != waves.end(); ++ptr, ++waveId) {

		auto& wave = **ptr;
		if (wave.tasks.size() == 0) {
			continue;
		}

		// update statistics
		auto& kernels = kernelSets[wave.category];
		tasksCount[wave.category] += wave.size();

		TIMER_MOVEON(denseTimer);

		// prepare stuff for dense calculation
		clCall(cl->mainDevice->mainQueue->enqueueWriteBuffer(buffer_tasks, CL_FALSE, 0, wave.size() * sizeof(PosteriorTask),  wave.tasks.data())); // start transfer
		PosteriorSchedule scheduling;
		scheduling.taskCount = wave.size();

		int groupSize = 0;
		if (wave.category == 0) {
			groupSize = mathex::ceilround(wave.maxHorizontal() + 1, 32);
		} else {
			float aux = (float)(wave.maxHorizontal() + 1) / this->maxWorkgroupSizes[wave.category];
			aux = (float)(wave.maxHorizontal() + 1) / std::ceil(aux);
			groupSize = mathex::ceilround((int)aux, 32);
		}

		int localMemory = PosteriorTask::localElements(groupSize, config->optimisation.useDoublePartition) * sizeof(buffer_type);

		TIMER_STOP(denseTimer);
		LOG_DEBUG << (wave.category == 0 ? "SHORT " : (wave.category == 1 ? "LONG  " : "LONGEST "))
			<< waveId << "/" << waves.size() << ": "
			<< "alloc-mb= " << MemoryTools::processCurrentVirtual() / 1000000
			<< ", sparse-mb= " << (sparseBytes[0] + sparseBytes[1]) / 1000000
			<< ", tasks= " << wave.tasks.size()
			<< ", group-size= " << groupSize
			<< ", max-h= " << wave.maxHorizontal()
			<< ", max-v= " << wave.maxVertical()
			<< ", local-mem= " << localMemory
			<< ", ";
		TIMER_MOVEON(denseTimer);

		// calculate dense posterior matrices
		LOG_DEBUG << "Pxy...";
		cl::NDRange globalRange(groupSize * scheduling.taskCount);
		cl::NDRange localRange(groupSize);
		cl::NDRange offsetRange(0);
		kernels(KernelSet::PARTITION_FORWARD).setArg(7, scheduling);
		kernels(KernelSet::PARTITION_FORWARD).setArg(8, cl::__local(sizeof(buffer_type) * kernels.localMemory(KernelSet::PARTITION_FORWARD, groupSize)));
		kernels(KernelSet::PARTITION_REVERSE).setArg(7, scheduling);
		kernels(KernelSet::PARTITION_REVERSE).setArg(8, cl::__local(sizeof(buffer_type) * kernels.localMemory(KernelSet::PARTITION_REVERSE, groupSize)));

		if (config->optimisation.divideHmmKernels) {
			kernels(KernelSet::HMM_FORWARD).setArg(7, scheduling);
			kernels(KernelSet::HMM_FORWARD).setArg(8,  cl::__local(sizeof(buffer_type) * kernels.localMemory(KernelSet::HMM_FORWARD, groupSize)));
			kernels(KernelSet::HMM_BACKWARD).setArg(7, scheduling);
			kernels(KernelSet::HMM_BACKWARD).setArg(8,  cl::__local(sizeof(buffer_type) * kernels.localMemory(KernelSet::HMM_BACKWARD, groupSize)));
			kernels(KernelSet::HMM_COMBINE).setArg(4, scheduling);
		} else {
			kernels(KernelSet::HMM_ALL).setArg(7, scheduling);
			kernels(KernelSet::HMM_ALL).setArg(8, cl::__local(sizeof(buffer_type) * kernels.localMemory(KernelSet::HMM_ALL, groupSize)));
		}

		kernels(KernelSet::FINALIZATION).setArg(4, scheduling);
		kernels(KernelSet::FINALIZATION).setArg(5, cl::__local(sizeof(buffer_type) * kernels.localMemory(KernelSet::FINALIZATION, groupSize)));

		cl::WaitForEvents(finishedSparse);

		clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
			kernels(KernelSet::PARTITION_FORWARD), offsetRange, globalRange, localRange, nullptr, phaseFinished[0])); // start after sparse kernels finish
		clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
			kernels(KernelSet::PARTITION_REVERSE), offsetRange, globalRange, localRange, nullptr, phaseFinished[1]));
		if (config->optimisation.divideHmmKernels) {
			clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
				kernels(KernelSet::HMM_FORWARD), offsetRange, globalRange, localRange, nullptr, phaseFinished[2]));
			clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
				kernels(KernelSet::HMM_BACKWARD), offsetRange, globalRange, localRange, nullptr, phaseFinished[3]));
			clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
				kernels(KernelSet::HMM_COMBINE), offsetRange, globalRange, localRange, nullptr, phaseFinished[4]));
		} else {
			clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
				kernels(KernelSet::HMM_ALL), offsetRange, globalRange, localRange, nullptr, phaseFinished[4]));
		}
		// wait for partition kernel to finish
		clCall(cl->mainDevice->mainQueue->enqueueNDRangeKernel(
			kernels(KernelSet::FINALIZATION), offsetRange, globalRange, localRange, nullptr, phaseFinished[5]));
		clCall(cl->mainDevice->mainQueue->enqueueReadBuffer(buffer_tasks, CL_TRUE, 0, wave.size() * sizeof(PosteriorTask), wave.tasks.data()));

		if (config->optimisation.kernelProfiling) {
			for (int i = 0; i < phaseFinished.size(); ++i) {
				timePhases[i] += clex::OpenCL::profileTimeMsec(*phaseFinished[i]);
			}
		}

		TIMER_STOP(denseTimer);
		LOG_DEBUG << "ok";

		// calculate sparse matrices
		TIMER_MOVEON(sparseTimer);

		// update sparse matrix offsets and write updated tasks
		size_t openClSparseElems = 0;
		for (auto& t : wave.tasks) {
			t.offset_Aij = openClSparseElems;
			size_t elemsSparse = OpenCLSparseMatrix::bytesNeeded(cms.lengths[t.i], t.numCells) / sizeof(buffer_type);
			size_t elemsSparseTransposed = OpenCLSparseMatrix::bytesNeeded(cms.lengths[t.j], t.numCells) / sizeof(buffer_type);
			openClSparseElems += std::max(elemsSparse, elemsSparseTransposed) + 4;
		}
		bool fitInHalf = 2 * openClSparseElems < auxiliaryElementsCount;

		LOG_DEBUG << ", fit=" << fitInHalf;

		clCall(cl->mainDevice->mainQueue->enqueueWriteBuffer(
			buffer_tasksSparse, CL_FALSE, 0, wave.size() * sizeof(PosteriorTask),  wave.tasks.data(), NULL, &finishedDense[0]));

		// execute kernels for both sparse and sparse transposed computation
		for (int qid = 0; qid < 2; ++qid) {

			// wait for unpack from previous wave to finish
			if (unpackThreads[qid].joinable()) { unpackThreads[qid].join(); }

			// this is important only when two sparse matrices do not fit in the auxiliary buffer
			if (qid == 1 && !fitInHalf && unpackThreads[0].joinable()) {
				LOG_DEBUG << "finishing...";
				unpackThreads[0].join();
			}

			string label = (qid == 0) ? ", Sxy..." : ", Syx...";
			int matrixHeight = (qid == 0)
				? (wave.category == 0 ? wave.maxVertical() + 1 : 4)
				: std::min(wave.maxHorizontal() + 1, this->maxWorkgroupSizes[wave.category]);

			LOG_DEBUG << label;
			localMemory = matrixHeight * sizeof(SparseMatrixType::index_type);
			cl::Buffer& buffer_currentAux = (qid == 1 && fitInHalf) ? subbuffer_auxiliaryHi : buffer_auxiliary;
			buffer_type* currentAuxiliary = auxiliaryMatrices + ((qid == 1 && fitInHalf) ? (auxiliaryElementsCount / 2) : 0);

			kernels(KernelSet::SPARSE).setArg(1, buffer_currentAux);
			kernels(KernelSet::SPARSE).setArg(4, scheduling);
			kernels(KernelSet::SPARSE).setArg(5, cl::__local(localMemory));
			kernels(KernelSet::SPARSE).setArg(6, qid);
			clCall(cl->mainDevice->queues[qid]->enqueueNDRangeKernel(kernels(KernelSet::SPARSE), offsetRange, globalRange, localRange, &finishedDense));
			clCall(cl->mainDevice->queues[qid]->enqueueReadBuffer(
				buffer_currentAux, CL_FALSE, 0, openClSparseElems * sizeof(buffer_type), currentAuxiliary, NULL, &finishedSparse[qid]));
			LOG_DEBUG << "sent...";

			// start new thread which generates and unpacks sparse matrices
			unpackThreads[qid] = std::thread(
				[qid, fillRate, this, currentAuxiliary, &finishedSparse, &wave, &matrices, &distances, &unpackTimer, &cms, &sparseBytes](){
				finishedSparse[qid].wait();
				//	TIMER_MOVEON(unpackTimer);
				int taskId;
				size_t bytes = 0;
				#pragma omp parallel for private(taskId) default(shared) schedule(dynamic) reduction(+:bytes)
				for(taskId = 0; taskId < wave.tasks.size(); taskId++) {
					auto& task = wave.tasks[taskId];
					int i = (qid == 0) ? task.i : task.j;
					int j = (qid == 0) ? task.j : task.i;
					SparseMatrixType* matrix = new SparseMatrixType(cms.lengths[i], cms.lengths[j]);
					matrices[i][j] = matrix;
					distances[i][j] = matrix->filterAndUnpack<OpenCLSparseMatrix>(
						currentAuxiliary + task.offset_Aij, static_cast<float>(fillRate));

					bytes += matrix->bytesNeeded();
				}
				sparseBytes[wave.category].fetch_add(bytes);
				//	TIMER_STOP(unpackTimer);
			});
		}

		LOG_DEBUG << endl;
	}

	LOG_DEBUG << endl;

	// join both unpacking threads
	if (unpackThreads[0].joinable()) { unpackThreads[0].join(); }
	if (unpackThreads[1].joinable()) { unpackThreads[1].join(); }
	if (cpuWaveThread.joinable()) { cpuWaveThread.join(); }

	if (config->optimisation.kernelProfiling) {
		std::for_each(phaseFinished.begin(), phaseFinished.end(), [](cl::Event* e) { delete e; });
	}

	delete [] auxiliaryMatrices;

	this->sparseBytes = sparseBytes[0] + sparseBytes[1] + sparseBytes[2];

	STATS_WRITE("memory.posterior MB", this->denseBytes / 1e6);
	STATS_WRITE("memory.sparse MB", this->sparseBytes/ 1e6);
	STATS_WRITE("memory.max posterior wave elems", waveMaxMemory->auxiliaryElemsCount);

	STATS_WRITE("scheduling.compact posterior tasks", tasksCount[0]);
	STATS_WRITE("scheduling.compact posterior waves", wavesCount[0]);
	STATS_WRITE("scheduling.long posterior tasks", tasksCount[1]);
	STATS_WRITE("scheduling.long posterior waves", wavesCount[1]);
	STATS_WRITE("scheduling.very long posterior tasks", tasksCount[2]);

	STATS_WRITE("scheduling.posterior cutoff", config->algorithm.posteriorCutoff);

	STATS_WRITE("time.1-1-task preparation", TIMER_SECONDS(prepareTimer));

	STATS_WRITE("time.1-2-a-partition forward", timePhases[0]);
	STATS_WRITE("time.1-2-b-partition reverse", timePhases[1]);
	STATS_WRITE("time.1-2-c-probabilistic forward", config->optimisation.divideHmmKernels ? timePhases[2] : 0);
	STATS_WRITE("time.1-2-d-probabilistic backward", config->optimisation.divideHmmKernels ? timePhases[3] : 0);
	STATS_WRITE("time.1-2-e-probabilistic combine", timePhases[4]);
	STATS_WRITE("time.1-2-f-finalisation", timePhases[5]);
	STATS_WRITE("time.1-2-gpu dense calculation", TIMER_SECONDS(denseTimer));

	STATS_WRITE("time.1-3-gpu sparse calculation", TIMER_SECONDS(sparseTimer));
	STATS_WRITE("time.1-4-sparse unpacking", TIMER_SECONDS(unpackTimer));
	STATS_WRITE("time.1-5-cpu processing", timeCpuProcessing);
}


void QuickPosteriorStage::computeWaveCpu(
	const PosteriorTasksWave& wave,
	const quickprobs::ContiguousMultiSequence& sequences,
	Array<float>& distances,
	Array<SparseMatrixType*>& matrices,
	std::atomic< ::size_t>& sparseBytes)
{
	if (wave.size() == 0) {
		return;
	}

	TIMER_CREATE(timer);
	TIMER_START(timer);

	// do all pairwise alignments for posterior probability matrices
	int numThreads = config->hardware.numThreads;
	int maxLayerSize = (wave.maxLength() + 1) * (wave.maxLength() + 1);
	std::vector<std::shared_ptr<BufferSet>> bufferSets(numThreads);
	for (auto& set : bufferSets) {
		set = std::shared_ptr<BufferSet>(new BufferSet(maxLayerSize));
	}

	int taskId;
	::size_t bytes = 0;
	#pragma omp parallel for private(taskId) reduction (+ : bytes) default(shared) schedule(dynamic)
	for(taskId = 0; taskId < wave.tasks.size(); taskId++) {
		auto& task = wave.tasks[taskId];

		// to make results identical as for cpu calculations restore original i,j coordinates
		int i = task.i;//(task.i < task.j) ? task.i : task.j;
		int j = task.j;//(task.i < task.j) ? task.j : task.i;

		int tid = omp_get_thread_num();
		BufferSet& set = *bufferSets[tid];

		computePairwise(*sequences.GetSequence(i), *sequences.GetSequence(j), set, distances[i][j]);

		auto m = new SparseMatrixType(sequences.lengths[i], sequences.lengths[j], set.f0(), config->algorithm.posteriorCutoff);
		auto mt = m->computeTranspose();

		bytes += m->bytesNeeded() + mt->bytesNeeded();

		matrices[i][j] = m;
		matrices[j][i] = mt;
		distances[j][i] = distances[i][j];
	}

	TIMER_STOP(timer);

	this->timeCpuProcessing += timer.seconds();
	sparseBytes = bytes;
}

