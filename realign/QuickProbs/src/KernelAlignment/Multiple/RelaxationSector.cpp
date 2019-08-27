#include "Common/MemoryTools.h"
#include "Common/Log.h"
#include "Common/mathex.h"

#include "Alignment/DataStructures/ContiguousMultiSequence.h"
#include "Alignment/DataStructures/PackedSparseMatrix.h"
#include "Alignment/DataStructures/SparseHelper.h"
#include "RelaxationSector.h"

#include "Common/deterministic_random.h"

using namespace quickprobs;

std::vector<std::shared_ptr<RelaxationSector>> RelaxationSector::generate(
	const Array<SparseMatrixType*>& sparseMatrices,
	const std::vector<unsigned int>& sortingMap,
	const Array<float>& distances,
	const Selectivity& selectivity,
	std::vector<unsigned int>& sparseOffsets,
	int sectorResolution)
{
	int numSeqs = sparseMatrices.size();
	int sectorSize =  mathex::ceildiv(numSeqs, sectorResolution);
	sectorResolution = mathex::ceildiv(numSeqs, sectorSize);

	sparseOffsets.resize(numSeqs * numSeqs);
	std::vector<std::shared_ptr<RelaxationSector>> sectors(sectorResolution * sectorResolution);
	
	std::vector<int> seeds(numSeqs * numSeqs);
	std::mt19937 eng;
	det_uniform_int_distribution<int> dist(0, RND_MAX);
	std::generate(seeds.begin(), seeds.end(), [&eng, &dist]()->int {
		return dist(eng);
	});

	// iterate through all sectors in parallel
	int sectorId;
	#pragma omp parallel for private(sectorId) default(shared) schedule(dynamic)
	for (sectorId = 0; sectorId < sectors.size(); ++sectorId) {

		int sector_i = sectorId / sectorResolution;
		int sector_j = sectorId % sectorResolution;

		int i_begin = sector_i * sectorSize;
		int j_begin = sector_j * sectorSize;
		int i_end = std::min(i_begin + sectorSize, numSeqs);
		int j_end = std::min(j_begin + sectorSize, numSeqs);
		int height = i_end - i_begin;
		int width = j_end - j_begin; 
			
		// number of tasks in sector depends on sector position:
		// - diagonal sector: only upper triangle,
		// - upper triangle sector: all tasks,
		// - bottom triangle sector: none tasks
		int taskCount = 0;
		if (sector_i == sector_j) {
			taskCount = width * (height - 1) / 2; // here width and height are equal to sectorSize
		} else if (sector_i < sector_j) {
			taskCount = width * height;
		}

		auto sector = std::shared_ptr<RelaxationSector>(new RelaxationSector(taskCount, sector_i, sector_j, i_begin, j_begin, height, width));
		sectors[sector_i * sectorResolution + sector_j] = sector;

	//	sumW = 0;

		// generate all possible tasks, ui and uj indicate unsorted indices
		auto task = sector->tasks.begin();
		
		for (int ui = i_begin; ui < i_end; ++ui) {
			for (int uj = j_begin; uj < j_end; ++uj) {
				int i = sortingMap[ui];
				int j = sortingMap[uj];

				if (ui < uj) {
					task->i = i;
					task->j = j;
					task->seed = seeds[i * numSeqs + j];
					task->acceptedCount = 0;

					int seed = task->seed;
					// check selectivity criterion
					for (int k = 0; k < distances.size(); ++k) {
						if (k != i && k != j) {	
							float d = selectivity.function(distances[i][k], distances[j][k]);
							seed = parkmiller(seed); // get next random number
							float v = selectivity.filter(selectivity.filter_a, selectivity.filter_b, d);
							float w = ((float)seed) * RND_MAX_INV - v;  

							if  (w < 0) {				
								++task->acceptedCount;
							} 
						}
					}

					++task;
				}
					
				if (i != j) {
					auto Sij = sparseMatrices[i][j];
					sparseOffsets[i * numSeqs + j] = sector->elementsCount;
					sector->elementsCount += mathex::ceilround(
						OpenCLSparseMatrix::bytesNeeded(Sij->getSeq1Length(), Sij->getNumCells()) / sizeof(buffer_type), 
						(::size_t)4);

					// update statistics
					sector->verticalLengths[Sij->getSeq1Length()] += 1;
					sector->horizontalLengths[Sij->getSeq2Length()] += 1;
					sector->sparseRows[Sij->getMaxRow()] += 1;
				} 
			}
		}
	}

	return sectors;
}

void quickprobs::RelaxationSector::packExternal(
	const Array<SparseMatrixType*>& sparseMatrices, 
	const std::vector<unsigned int>& sparseOffsets,
	const std::vector<unsigned int>& sortingMap,
	buffer_type* rawData)
{	
	//LOG_DEBUG << "packing (" << this->sector_i << "," << sector_j << ")...";
	if (this->rawData == nullptr) {
	//	LOG_DEBUG << (int64_t)rawData << " finished!" << std::endl;
		this->rawData = rawData;
		pack(sparseMatrices, sparseOffsets, sortingMap, rawData);
	} else {
	//	LOG_DEBUG << "already packed!" << std::endl;
	}
}

void quickprobs::RelaxationSector::packInternal(
	const Array<SparseMatrixType*>& sparseMatrices, 
	const std::vector<unsigned int>& sparseOffsets, 
	const std::vector<unsigned int>& sortingMap)
{
	if (!internallyAllocated) { // check if internal buffer already set
		rawData = new buffer_type[elementsCount];
		internallyAllocated = true;
		this->pack(sparseMatrices, sparseOffsets, sortingMap, rawData);
	} 
}


void quickprobs::RelaxationSector::pack(
	const Array<SparseMatrixType*>& sparseMatrices, 
	const std::vector<unsigned int>& sparseOffsets,
	const std::vector<unsigned int>& sortingMap,
	buffer_type* rawData)
{	
	// pack only when necessary
	int numSeqs = sparseMatrices.size();
	int matrixId;

#pragma omp parallel for private(matrixId) default(shared) schedule(dynamic)
	for (matrixId = 0; matrixId < width * height; ++matrixId) {

		int ui = i_begin + matrixId / width;
		int uj = j_begin + matrixId % width;

		int i = sortingMap[ui];
		int j = sortingMap[uj];

		if (i != j) {
			auto Sij = sparseMatrices[i][j];
			Sij->pack<OpenCLSparseMatrix>(rawData + sparseOffsets[i * numSeqs + j], false, 1.0f);
		}
	}
}


void quickprobs::RelaxationSector::unpack(
	const buffer_type* buffer, 
	const std::vector<unsigned int>& sparseOffsets, 
	Array<SparseMatrixType*>& sparseMatrices) const
{
	int numSeqs = sparseMatrices.size();
	int taskId;

	#pragma omp parallel for private(taskId) default(shared) schedule(dynamic)
	for(taskId = 0; taskId < tasks.size(); taskId++) {
		auto& task = tasks[taskId];
		int i = task.i;
		int j = task.j;

		auto matrix = sparseMatrices[i][j];
		matrix->unpack<OpenCLSparseMatrix>(buffer + sparseOffsets[i * numSeqs + j]);
		sparseMatrices[j][i]->fillTransposed(*matrix);
	}
}


