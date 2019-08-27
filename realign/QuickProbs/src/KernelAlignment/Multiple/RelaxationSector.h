#pragma once
#include <vector>
#include <list>
#include <map>
#include <functional>

#include "Alignment/Multiple/Selectivity.h"
#include "Alignment/DataStructures/SparseMatrixType.h"
#include "Alignment/Multiple/AuxiliaryStructures.h"
#include "ITasksWave.h"

#include "Common/Array.h"

namespace quickprobs 
{
class ContiguousMultiSequence;

class RelaxationSector : public ITasksWave
{
public:
	
	const std::vector<RelaxationTask>& getTasks() const { return tasks; }

	int minSparseRow() const { return sparseRows.cbegin()->first; }
	int maxSparseRow() const { return sparseRows.crbegin()->first; }
	buffer_type* getRawData() { return rawData; }
	::size_t getElementsCount() { return elementsCount; }
	::size_t getBytes() { return elementsCount * sizeof(buffer_type); }

	virtual ::size_t size() const { return tasks.size(); }

	void resetRawData() { 
	//	LOG_DEBUG << "reset (" << sector_i << "," << sector_j << ")" << std::endl;
		if (internallyAllocated) { delete [] rawData; }
		rawData = nullptr; 
	}

	static std::vector<std::shared_ptr<RelaxationSector>> generate(
		const Array<SparseMatrixType*>& sparseMatrices,
		const std::vector<unsigned int>& sortingMap,
		const Array<float>& distances,
		const Selectivity& selectivity,
		std::vector<unsigned int>& sparseOffsets,

		int sectorResolution);

	RelaxationSector(int taskCount, int sector_i, int sector_j, int i_begin, int j_begin, int height, int width) 
		: ITasksWave(), tasks(taskCount), sector_i(sector_i), sector_j(sector_j), 
		i_begin(i_begin), j_begin(j_begin), height(height), width(width), 
		elementsCount(0), rawData(nullptr), internallyAllocated(false)
	{}
	
	~RelaxationSector() 
	{ 
		if (internallyAllocated) { delete [] rawData; }
	}
	
	void packExternal(
		const Array<SparseMatrixType*>& sparseMatrices, 
		const std::vector<unsigned int>& sparseOffsets,
		const std::vector<unsigned int>& sortingMap,
		buffer_type* rawData);

	void packInternal(
		const Array<SparseMatrixType*>& sparseMatrices, 
		const std::vector<unsigned int>& sparseOffsets,
		const std::vector<unsigned int>& sortingMap);
	
	void unpack(
		const buffer_type* buffer,
		const std::vector<unsigned int>& sparseOffsets,
		Array<SparseMatrixType*>& sparseMatrices) const;

protected:
	int sector_i;
	int sector_j;
	int i_begin;
	int j_begin;
	int height;
	int width;
	
	std::vector<RelaxationTask> tasks;
	buffer_type* rawData;
	bool internallyAllocated;
	std::map<int, int> sparseRows;
	::size_t elementsCount;

	void pack(
		const Array<SparseMatrixType*>& sparseMatrices, 
		const std::vector<unsigned int>& sparseOffsets,
		const std::vector<unsigned int>& sortingMap,
		buffer_type* rawData);
};

};
