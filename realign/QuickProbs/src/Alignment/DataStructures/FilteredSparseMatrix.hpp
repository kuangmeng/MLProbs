#pragma once
#include "FilteredSparseMatrix.h"


template <class CellType>
float FilteredSparseMatrix::calculateCutoff(const CellType* data, int count, int toRetain, int& newCount) const
{
	// iterate over matrices
	std::vector<int> buckets(cutoffs.size() + 1);
	int cutoffIndex = cutoffs.size() - 1;

	if (count > toRetain) {
		std::fill(buckets.begin(), buckets.end(), 0);
		for (int e = 0; e < count; ++e) {
			int c = cutoffs.size() / 2;
			int h = c / 2;
			while (h) {
				if (data[e].getValue() < cutoffs[c]) {	c += h; } 
				else { c -= h; }
				h >>= 1;
			}

			++buckets[c];		
		}

		// iterate over buckets to determine cutoff
		newCount = 0;

		for (int i = 0; i < buckets.size(); ++i) {
			newCount += buckets[i];
			if (newCount > toRetain) {
				newCount -= buckets[i];
				cutoffIndex = std::max(i - 1, 0);
				break;
			}
		}
	}

	if (newCount == 0) {
		newCount = count;
		return cutoffs.back();
	}

	return cutoffs[cutoffIndex];
}


template <class InMatrixType>
float FilteredSparseMatrix::filterAndUnpack(const buffer_type* ptr, float fillRate)
{
	const float* temp = (float*)ptr;
	float distance = temp[0];  
	int count = (int)temp[1];
	
	::size_t minElements = std::min(width, height);
	::size_t extraElements = width * height - minElements;
	::size_t toRetain = (::size_t)(extraElements * fillRate) + minElements;

	if (toRetain >= count) {
		return this->unpack<InMatrixType>(ptr);
	} else {
		const typename InMatrixType::index_type* rowSizesIn = InMatrixType::pointerToRowSizes(ptr) ;
		const typename InMatrixType::cell_type* dataIn = InMatrixType::pointerToData(ptr, getSeq1Length());

		float cutoff = calculateCutoff(dataIn, count, toRetain, numCells);
		::size_t allocSize = this->bytesNeeded();		
		assert(numCells >= 0 && numCells <= height * width);
	
		// allocate buffer if necessary
		if (buffer == nullptr) {
			buffer = new buffer_type[allocSize / sizeof(buffer_type)];
		} 

		// set pointers
		this->rowSizes = pointerToRowSizes(buffer);
		this->rowIndices = pointerToRowIndices(buffer, getSeq1Length());
		this->data = pointerToData(buffer, getSeq1Length());

		rowSizes[0] = -1;
		rowIndices[0] = 0;

		const typename InMatrixType::cell_type* currentIn = dataIn;
		cell_type* currentOut = data;

		for (int i = 1; i < height; ++i) {
			int newSize = 0;
			rowIndices[i] = currentOut - this->data;
		
			for (int j = 0; j < rowSizesIn[i]; ++j) {
				if (currentIn->getValue() >= cutoff) {
					currentOut->setColumn(currentIn->getColumn());
					currentOut->setValue(currentIn->getValue());
					++currentOut;
					++newSize;
				}
				++currentIn;
			}

			rowSizes[i] = newSize;
		}

		this->rowSizes[0] = *std::max_element(this->rowSizes + 1, this->rowSizes + height);
		int maxrow = getMaxRow();	
	
		if (maxrow > width) {
			std::cout << "Max row size error!" << std::endl;
			throw std::runtime_error("Max row size error");
		}
	
		assert(getNumCells() == std::accumulate(this->rowSizes + 1, this->rowSizes + height, 0));
		assert(maxrow < width);
		
		return distance;
	}
}