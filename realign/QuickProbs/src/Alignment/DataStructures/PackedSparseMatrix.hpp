#pragma once
#include "PackedSparseMatrix.h"

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstring>

template <class InMatrixType>
float PackedSparseMatrix::unpack(const buffer_type* ptr)
{
	const float* temp = (float*)ptr;
	
	float distance = temp[0];  
	numCells = (int)temp[1];
	
	assert(numCells >= 0 && numCells <= height * width);
	
	// allocate buffer if necessary
	::size_t allocSize = this->bytesNeeded();
	if (buffer == nullptr) {
		buffer = new buffer_type[allocSize / sizeof(buffer_type)];
	} 


	// copy data and set pointers
	this->rowSizes = pointerToRowSizes(buffer);
	this->rowIndices = pointerToRowIndices(buffer, getSeq1Length());
	this->data = pointerToData(buffer, getSeq1Length());
	
	const typename InMatrixType::cell_type* inData = InMatrixType::pointerToData(ptr, getSeq1Length());

	std::memcpy(buffer, ptr, this->metaBytesNeeded());
	std::transform(inData, inData + numCells, this->data, [](typename InMatrixType::cell_type in)->cell_type{
		cell_type out;
		out.setColumn(in.getColumn());
		out.setValue(in.getValue());
		return out;
	});

	this->rowSizes[0] = *std::max_element(this->rowSizes + 1, this->rowSizes + height);
	
	if (getMaxRow() > width) {
		std::cout << "Max row size error: " << getMaxRow() << std::endl;
		throw std::runtime_error("Max row size error");
	}
	
	assert(getNumCells() == std::accumulate(this->rowSizes + 1, this->rowSizes + height, 0));
	assert(getMaxRow() < width);

	return distance;
}

template <class OutMatrixType>
void PackedSparseMatrix::pack(buffer_type *ptr, bool freeBuffer, float initialDistance)
{
	// copy metadata
	std::memcpy(ptr, this->buffer, this->metaBytesNeeded());

	// copy sparse elements
	typename OutMatrixType::cell_type* outData = OutMatrixType::pointerToData(ptr, getSeq1Length());

	std::transform(this->data, this->data + numCells, outData, [](const cell_type& in)->typename OutMatrixType::cell_type{
		typename OutMatrixType::cell_type out;
		out.setColumn(in.getColumn());
		out.setValue(in.getValue());
		return out;
	});

	((float*)ptr)[0] = initialDistance;
	((float*)ptr)[1] = getNumCells(); 

	// effectively clear the memory
	if (freeBuffer) {
		delete [] buffer;
		buffer = nullptr;
	}
}