#pragma once
#include <vector>
#include <list>
#include <map>
#include <memory>

#include "ITasksWave.h"
#include "Alignment/Multiple/AuxiliaryStructures.h"
#include "Common/Array.h"

namespace quickprobs
{
class ContiguousMultiSequence;

class PosteriorTasksWave : public ITasksWave
{
public:
	std::vector<PosteriorTask> tasks;
	::size_t posteriorElemsCount;
	::size_t auxiliaryElemsCount;
	::size_t metaBytes;
	::size_t minSparseBytes;

	int category;
	
	virtual ::size_t size() const { return tasks.size(); } 

	PosteriorTasksWave(int category) 
		: posteriorElemsCount(0), auxiliaryElemsCount(0), minSparseBytes(0), metaBytes(0), category(category) {}

	static std::map<int, std::vector<PosteriorTask>> generateTasks(
		const ContiguousMultiSequence& sequences,
		const Array<float>& pids,
		::size_t maxWaveBytes, 
		::size_t maxAllocationBytes, 
		int lengthThreshold);

	static std::list<std::shared_ptr<PosteriorTasksWave>> generate(
		std::vector<PosteriorTask>& tasks,
		::size_t maxWaveBytes,		
		::size_t maxAllocationBytes,
		int columnElements,
		int category);

};

};