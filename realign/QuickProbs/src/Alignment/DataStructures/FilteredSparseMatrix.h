#pragma once 
#include "PackedSparseMatrix.h"

class FilteredSparseMatrix : public PackedSparseMatrix {
public:
	
	static void initialise(float initialCutoff, int binsCount);

	FilteredSparseMatrix(int seq1Length, int seq2Length)
		: PackedSparseMatrix(seq1Length, seq2Length) {}
	
	FilteredSparseMatrix(int seq1Length, int seq2Length, ::size_t numCells) 
		: PackedSparseMatrix(seq1Length, seq2Length, numCells) {} 

	FilteredSparseMatrix(int seq1Length, int seq2Length, const float*  posterior, float cutoff)
		: PackedSparseMatrix(seq1Length, seq2Length, posterior, cutoff) {}

	FilteredSparseMatrix* computeTranspose() const;

	template <class InMatrixType>
	float filterAndUnpack(const buffer_type* ptr, float fillRate);

	template <class CellType>
	float calculateCutoff(const CellType* data, int count, int toRetain, int& newCount) const;

	int filter(float cutoff);

	int filterLongRows();

protected:
	static std::vector<float> cutoffs;

};

#include "FilteredSparseMatrix.hpp"