#include <functional>
#include <algorithm>

#include "Alignment/DataStructures/JaggedMatrix.h"
#include "Alignment/DataStructures/Sequence.h"
#include "Alignment/DataStructures/MultiSequence.h"
#include "Alignment/DataStructures/ContiguousMultiSequence.h"
#include "Alignment/DataStructures/SparseMatrixType.h"
#include "PosteriorTasksWave.h"

using namespace quickprobs;


std::map<int, std::vector<PosteriorTask>> quickprobs::PosteriorTasksWave::generateTasks(
	const ContiguousMultiSequence& sequences, 
	const Array<float>& pids, 
	::size_t maxWaveBytes, 
	::size_t maxAllocationBytes, 
	int lengthThreshold)
{
	int numSeq = sequences.count();

	// at first generate all possible tasks and compute corresponding layer sizes
	std::map<int, std::vector<PosteriorTask>> out;
	out[0].resize(numSeq * (numSeq - 1) / 2);
	out[1].resize(numSeq * (numSeq - 1) / 2);
	out[2].resize(numSeq * (numSeq - 1) / 2);
	
	std::vector<std::vector<PosteriorTask>::iterator> currentTasks(3);
	currentTasks[0] = out[0].begin();
	currentTasks[1] = out[1].begin();
	currentTasks[2] = out[2].begin();

	for (int i = 0; i < numSeq; ++i) {
		for (int j = i + 1; j < numSeq; ++j) {
			int vertical = sequences.lengths[i] >= sequences.lengths[j] ? i : j;
			int horizontal = sequences.lengths[i] >= sequences.lengths[j] ? j : i;
			
			int height = sequences.lengths[vertical] + 1;
			int width = sequences.lengths[horizontal] + 1;
			
			std::vector<PosteriorTask>::iterator* task;
			
			int layerBytes = JaggedMatrix<float>::size(width, height, jagSize_float) * sizeof(buffer_type);

			// if too large for GPU execution
			if (4 * layerBytes > maxWaveBytes || 3 * layerBytes > maxAllocationBytes) {
				task = &currentTasks[2]; // very long task
			} else if (width < lengthThreshold) {
				task = &currentTasks[0]; // short task
			} else {
				task = &currentTasks[1]; // long task
			}

			(**task).i = vertical;
			(**task).j = horizontal;
			(**task).pid = pids[i][j];
			(**task).width = width;
			(**task).height = height;
			(**task).layerSize = JaggedMatrix<float>::size(width, height, jagSize_float);
			++(*task);
		}
	}

	// erase empty slots
	out[0].erase(currentTasks[0], out[0].end());
	out[1].erase(currentTasks[1], out[1].end());
	out[2].erase(currentTasks[2], out[2].end());

	return out;
}


std::list<std::shared_ptr<PosteriorTasksWave>> quickprobs::PosteriorTasksWave::generate(
	std::vector<PosteriorTask>& tasks, 
	::size_t maxWaveBytes, 
	::size_t maxAllocationBytes,
	int columnElements,
	int category)
{
	std::list<std::shared_ptr<PosteriorTasksWave>> waves;
	
	if (tasks.size() == 0) {
		return waves;
	}

	// sort tasks by widths sizes
	std::stable_sort(tasks.begin(), tasks.end(), [](const PosteriorTask &t1, const PosteriorTask &t2) { 
		return t1.width > t2.width; 
	});

	// generate set of waves
	::size_t freeElemsCount = maxWaveBytes / sizeof(buffer_type);
	::size_t freeAllocationElems = maxAllocationBytes / sizeof(buffer_type);
	auto low = tasks.begin();
	std::shared_ptr<PosteriorTasksWave> wave(new PosteriorTasksWave(category));

	for (auto task = tasks.begin(); task != tasks.end(); ++task) {
		if (task != low && task->height < wave->maxHorizontal()) { // this is correct as widest task is the first one
			// transpose matrix and recalculate
			std::swap(task->i, task->j);
			std::swap(task->width, task->height);
			task->layerSize = JaggedMatrix<float>::size(task->width, task->height, jagSize_float); // use height as width and opposite
		}

		const auto& len_i = task->height - 1;
		const auto& len_j = task->width - 1;
		
		// this is for both Sxy and Syx
		wave->minSparseBytes += 
			SparseMatrixType::BaseType::bytesNeeded(len_i, std::min(task->height, task->width)) + 
			SparseMatrixType::BaseType::bytesNeeded(len_j, std::min(task->height, task->width));

		// reserve layers for auxiliary buffer and output sparse matrices
		int auxiliaryRequirement = std::max(
			(::size_t)task->layerSize * 3,
			SparseMatrixType::BaseType::bytesNeeded(len_i, len_i * len_j) / sizeof(buffer_type));
		int posteriorRequirement = task->layerSize;

		// make auxiliary buffer 8-byte aligned
		if (auxiliaryRequirement % 2 != 0) {
			auxiliaryRequirement++;
		}

		int totalRequirement = posteriorRequirement + auxiliaryRequirement + sizeof(PosteriorTask) + columnElements;  

		// if there is no place for next task
		if ((freeElemsCount < totalRequirement) || (freeAllocationElems < auxiliaryRequirement)) {
			wave->tasks.insert(wave->tasks.begin(), low, task); // copy tasks to existing wave
			waves.push_back(wave);

			// create new wave
			wave = std::shared_ptr<PosteriorTasksWave>(new PosteriorTasksWave(category)); // create new wave
			freeElemsCount = maxWaveBytes / sizeof(buffer_type);
			freeAllocationElems = maxAllocationBytes / sizeof(buffer_type);
			low = task;
		}

		// update length statistics if necessary
		wave->verticalLengths[len_i] += 1;
		wave->horizontalLengths[len_j] += 1;

		task->offset_Aij = wave->auxiliaryElemsCount; // set offsets in buffers
		task->offset_Pij = wave->posteriorElemsCount;

		wave->auxiliaryElemsCount += auxiliaryRequirement; // increment sizes of wave buffers
		wave->posteriorElemsCount += posteriorRequirement;

		freeElemsCount -= totalRequirement;
		freeAllocationElems -= auxiliaryRequirement; // largest buffer to be allocated
	}

	wave->tasks.insert(wave->tasks.begin(), low, tasks.end()); // copy last tasks to wave
	waves.push_back(wave);

	return waves;
}


