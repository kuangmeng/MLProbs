#include "PackedSparseMatrix.h"

#include <algorithm>
#include <numeric>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cstring>
#include <limits>


int PackedSparseMatrix::sparseThreshold = std::numeric_limits<int>::max();

float PackedSparseMatrix::getValue (int row, int col) const
{
	for (int i = 0; i < rowSizes[row]; i++) {
		if (data[rowIndices[row] + i].getColumn() == col) {
			return data[rowIndices[row] + i].getValue();
		}
	}
	return 0;
}

PackedSparseMatrix::PackedSparseMatrix(int seq1Length, int seq2Length, ::size_t numCells) :
	height(seq1Length + 1), width(seq2Length + 1), numCells(numCells), buffer(nullptr)
{
	// allocate memory
	//if (numCells > 0) {
		buffer = new buffer_type[BaseType::bytesNeeded(seq1Length, numCells) / sizeof(buffer_type)];

		this->rowSizes = pointerToRowSizes(buffer);
		this->rowIndices = pointerToRowIndices(buffer, seq1Length);
		this->data = pointerToData(buffer, seq1Length);
		rowSizes[0] = -1;
		rowIndices[0] = 0; 
	//}
}

PackedSparseMatrix::PackedSparseMatrix(int seq1Length, int seq2Length, const float* posterior, float cutoff) :
	height(seq1Length + 1), width(seq2Length + 1)
{
	int size = (seq1Length + 1) * (seq2Length + 1); 
	numCells = std::count_if(posterior, posterior + size, [&cutoff](float x)->bool {
		return x >= cutoff;
	});

	// allocate memory
	buffer = new buffer_type[bytesNeeded() / sizeof(buffer_type)];
	
	this->rowSizes = pointerToRowSizes(buffer);
	this->rowIndices =pointerToRowIndices(buffer, seq1Length);
	this->data = pointerToData(buffer, seq1Length);
	rowSizes[0] = -1;
	rowIndices[0] = 0;

	// build sparse matrix
	const float* postPtr = posterior + width;				// note that we're skipping the first row here

	auto dataPtr = data;
	for (int i = 1; i < height; i++){
		postPtr++;                                            // and skipping the first column of each row
		rowIndices[i] = dataPtr - data; // added

		int rowSpace = 0;

		for (int j = 1; j < width; j++) {
			if (*postPtr >= cutoff) {
				dataPtr[rowSpace].setColumn(j);
				dataPtr[rowSpace].setValue(*postPtr);
				rowSpace++;
			}
			postPtr++;
		}

		rowSizes[i] = rowSpace;
		dataPtr += rowSpace;
	}

	this->rowSizes[0] = *std::max_element(this->rowSizes + 1, this->rowSizes + height);
	assert(getMaxRow() < width);
	assert(getNumCells() == std::accumulate(this->rowSizes + 1, this->rowSizes + height, 0));
}

PackedSparseMatrix::~PackedSparseMatrix()
{
	if (buffer != nullptr) {
		delete [] buffer;
	}
}


PackedSparseMatrix* PackedSparseMatrix::computeTranspose() const
{
	// create a new sparse matrix
	PackedSparseMatrix *out = new PackedSparseMatrix(this->getSeq2Length(), this->getSeq1Length(), this->getNumCells());
	out->fillTransposed(*this);
	return out;
}

void PackedSparseMatrix::fillTransposed(const PackedSparseMatrix& other)
{
	this->numCells = other.numCells;
	rowSizes[0] = -1;
	rowIndices[0] = 0; 

	// compute row sizes
	std::fill(rowSizes, rowSizes + height, 0);
	//std::fill(rowIndices, rowIndices + height, 0);

	for (int i = 0; i < numCells; ++i) { ++rowSizes[other.data[i].getColumn()]; }

	rowIndices[0] = numCells; 
	rowIndices[1] = 0;

	// compute row indices 
	for (int i = 2; i < height; i++) {
		rowIndices[i] = rowIndices[i-1] + rowSizes[i-1];
	}

	// reuse sizes of storing offsets
	int* offsets = rowSizes; 
	std::fill(offsets, offsets + height, 0);
	
	// now fill in data
	for (int i = 1; i < other.height; ++i) {
		for (int j = 0; j < other.rowSizes[i]; ++j) {
			const auto& cell = other.data[other.rowIndices[i] + j];

			data[rowIndices[cell.getColumn()] + offsets[cell.getColumn()]].setColumn(i);
			data[rowIndices[cell.getColumn()] + offsets[cell.getColumn()]].setValue(cell.getValue());

			++offsets[cell.getColumn()];
		}
	}

	this->rowSizes[0] = *std::max_element(this->rowSizes + 1, this->rowSizes + height);
	assert(getMaxRow() < width);
	assert(getNumCells() == std::accumulate(this->rowSizes + 1, this->rowSizes + height, 0));
}


void PackedSparseMatrix::shrinkToFit()
{
	buffer_type* newBuffer = new buffer_type[this->bytesNeeded() / sizeof(buffer_type)];
	std::memcpy(newBuffer, buffer, this->bytesNeeded());
	delete [] buffer;
	buffer = newBuffer;
}


std::string PackedSparseMatrix::toString()
{
	std::ostringstream oss;
	oss << "numCells = " << numCells << std::endl;
	for (int i = 0; i < height; ++i) {
		oss << "sizes = " << rowSizes[i] << "(i), " << *((float*)&rowSizes[i]) << "(f), " 
			<< "indices = " << rowIndices[i] << "(i), " << *((float*)&rowIndices[i]) << "(f), " << std::endl;
	}

	for (int i = 0; i < numCells; ++i) {
		oss << "cells = " << "(" << data[i].getColumn() << ", " << data[i].getValue() << "), "; 
	}

	return oss.str();
}

std::vector<float> * PackedSparseMatrix::getPosterior() const
{
	// create a new posterior matrix
	std::vector<float> *posteriorPtr = new std::vector<float>(height * width, 0);
	float *row = posteriorPtr->data() + width;

	// build the posterior matrix
	for (int i = 1; i < height; i++, row += width){
		for (int j = 0; j < rowSizes[i]; j++){
			auto elem = data[rowIndices[i] + j];
			row[elem.getColumn()] = elem.getValue();
		}
	}

	return posteriorPtr;
}

