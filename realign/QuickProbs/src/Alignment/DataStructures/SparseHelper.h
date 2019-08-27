#pragma once
#include <map>
#include "Common/mathex.h"
#include "SparseMatrixType.h"

#include "Common/Array.h"

class SparseHelper {
public:
	
	static ::size_t totalElements(Array<SparseMatrixType*>& matrices);
	static double sumOfElements(Array<SparseMatrixType*>& matrices);
	static double sumOfElements(SparseMatrixType& matrix);

	static void getHistogram(SparseMatrixType& matrix, std::map<int,int>& histogram);

	static void getHistogram(Array<SparseMatrixType*>& matrices, std::map<int,int>& histogram);

	static void filterStatic(Array<SparseMatrixType*>& matrices, float cutoff);

	static void filterLocal(
		Array<SparseMatrixType*>& matrices, 
		float initialCutoff,
		float remainFraction, 
		int maxSparseRow, 
		bool sparseReference);

	static void filterLocal(Array<SparseMatrixType*>& matrices);
	
	static int filterGlobal(Array<SparseMatrixType*>& matrices, float remainFraction, bool sparseReference);
};