#pragma once
#include "SparseEntry.h"

#include <CL/cl.h>
#include <vector>
#include <cstdint>

typedef int32_t buffer_type;


template <class IndexType, class ColumnType, class ValueType>
class SparseMatrixBase
{
public:
	typedef IndexType index_type;
	typedef ColumnType column_type;
	typedef ValueType value_type;
	typedef SparseEntry<column_type, value_type> cell_type;

	static const ::size_t metaBytesNeeded(int seq1Length) {
		return 2 * sizeof(float) +	2 * (seq1Length + 1) * sizeof(index_type);	
	}

	static ::size_t bytesNeeded(int seq1Length, ::size_t numCells) {
		return metaBytesNeeded(seq1Length) + numCells * sizeof(cell_type);
	}

	static index_type* pointerToRowSizes(buffer_type* ptr) {
		return (index_type*)((float*)ptr + 2);
	}

	static const index_type* pointerToRowSizes(const buffer_type* ptr) {
		return (const index_type*)((float*)ptr + 2);
	}

	static index_type* pointerToRowIndices(buffer_type* ptr, int seq1Length) {
		return (index_type*)((float*)ptr + 2) + seq1Length + 1;
	}

	static const index_type* pointerToRowIndices(const buffer_type* ptr, int seq1Length) {
		return (const index_type*)((float*)ptr + 2) + seq1Length + 1;
	}

	static cell_type* pointerToData(buffer_type* ptr, int seq1Length) {
		return (cell_type*)(ptr + metaBytesNeeded(seq1Length) / sizeof(buffer_type));
	}

	static const cell_type* pointerToData(const buffer_type* ptr, int seq1Length) {
		return (const cell_type*)(ptr + metaBytesNeeded(seq1Length) / sizeof(buffer_type));
	}
};

class FixedSparseMatrix : public SparseMatrixBase<int, uint16_t, uint16_t> {};


//typedef SparseMatrixBase<int, int, float> FloatSparseMatrix;
//typedef SparseMatrixBase<int, uint16_t, uint16_t> FixedSparseMatrix;
typedef SparseMatrixBase<cl_int, cl_int, cl_float> OpenCLSparseMatrix;
