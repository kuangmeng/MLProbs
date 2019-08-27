#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>

#include "Common/mathex.h"
#include "SparseHelper.h"


::size_t SparseHelper::totalElements(Array<SparseMatrixType*>& matrices)
{
	::size_t elems = 0; 
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = i + 1; j < matrices.size(); ++j) {
			elems += matrices[i][j]->getNumCells();
		}
	}
	return elems;
}


double SparseHelper::sumOfElements(Array<SparseMatrixType*>& matrices)
{
	double sum = 0.0;
	::size_t elems = 0; 
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = i + 1; j < matrices.size(); ++j) {
			sum += sumOfElements(*matrices[i][j]);
		}
	}
	return sum;
}

double SparseHelper::sumOfElements(SparseMatrixType& matrix)
{
	double sum = 0.0;
	for (int i = 1; i < matrix.getHeight(); ++i) {
		for (int j = 0; j < matrix.getRowSize(i); ++j) {
			sum += matrix.getValue(i, j);
		}
	}
	return sum;
}



void SparseHelper::filterLocal( 
	Array<SparseMatrixType*>& matrices, 
	float initialCutoff,
	float remainFraction, 
	int maxSparseRow, 
	bool sparseReference )
{
	float multipiler = 1.1;
	int cutoffsCount = std::log(1.0 / initialCutoff) / std::log(multipiler) + 1;
	
	std::vector<float> cutoffs(cutoffsCount);
	std::vector<int> buckets(cutoffs.size() + 1);
	cutoffs[cutoffs.size() - 1] = initialCutoff;
	
	for (int i = cutoffs.size() - 1; i > 0; --i) {
		cutoffs[i-1] = cutoffs[i] * multipiler;
	}

	// iterate over matrices
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = 0; j < matrices.size(); ++j) {
			
			if (i == j) { continue; } // omit diagonal elements
			
			auto& matrix = *dynamic_cast<SparseMatrixType*>(matrices[i][j]);
			::size_t totalSparseElements = matrix.getNumCells();
			::size_t totalDenseElements = matrix.getWidth() * matrix.getHeight();
			::size_t totalElements = sparseReference ? totalSparseElements : totalDenseElements;
			::size_t toRetain = (::size_t)(totalElements * remainFraction);

			int cutoffIndex = 0;

			if (matrix.getNumCells() < toRetain) { 
				
				if (matrix.getMaxRow() < maxSparseRow) {
					// if no filtering is needed - move to another matrix
					continue;
				} else {
					// select one before last cutoff
					cutoffIndex = cutoffs.size() - 2;
				}
			} else {

				std::fill(buckets.begin(), buckets.end(), 0);
				auto matrixData = matrix.getData();
				for (int e = 0; e < matrix.getNumCells(); ++e) {
					int c;
					for (c = 0; c < cutoffs.size(); ++c) {
						if (matrixData[e].getValue() > cutoffs[c]) {
							break;
						}
					}

					++buckets[c];		
				}

				// iterate over buckets to determine cutoff
				int added = 0;
			
				for (int i = 0; i < buckets.size(); ++i) {
					added += buckets[i];
					if (added > toRetain) {
						cutoffIndex = i;
						break;
					}
				}
			}

			do {
		//		matrix.filter(cutoffs[cutoffIndex]);
				--cutoffIndex;
			} while (cutoffIndex >= 0 && matrix.getMaxRow() > maxSparseRow);

			//matrix.shrinkToFit();
		}
	}
}

void SparseHelper::filterLocal( Array<SparseMatrixType*>& matrices )
{
	// iterate over matrices
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = 0; j < matrices.size(); ++j) {

			if (i == j) { continue; } // omit diagonal elements

			auto& matrix = *dynamic_cast<SparseMatrixType*>(matrices[i][j]);

			//matrix.filter(); // filter using last stored cutoff
		}
	}
}

int SparseHelper::filterGlobal( Array<SparseMatrixType*>& matrices, float remainFraction, bool sparseReference )
{
	// generate histogram
	std::vector<float> cutoffs(30);
	std::vector<int> buckets(cutoffs.size() + 1);
	cutoffs[0] = 1.0;
	
	for (int i = 1; i < cutoffs.size(); ++i) {
		cutoffs[i] = cutoffs[i-1] / 1.2f;
	}

	::size_t totalSparseElements = 0;
	::size_t totalDenseElements = 0;

	// iterate over matrices
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = 0; j < matrices.size(); ++j) {
			
			if (i == j) {
				continue;
			}
			auto& matrix = *matrices[i][j];
			totalSparseElements += matrix.getNumCells();
			totalDenseElements += matrix.getHeight() * matrix.getWidth();

			for (auto elem = matrix.getRowPtr(0); elem < matrix.getRowPtr(0) + matrix.getNumCells() ; ++elem) {
				// iterate over borders
				bool bucketFound = false;
				for (int c = 0; c < cutoffs.size(); ++c) {
					
					if (elem->getValue() > cutoffs[c]) {
						++buckets[c];
						bucketFound = true;
						break;
					}
				}

				if (bucketFound == false) {
					++buckets[cutoffs.size()];
				}
			}
		}
	}

	::size_t totalElements = sparseReference ? totalSparseElements : totalDenseElements;

	int toRetain = static_cast<int>(totalElements * remainFraction); // compute number of elements to retain
	float cutoff = 0;

	// iterate over buckets to determine cutoff
	int current = 0;
	for (int i = 0; i < buckets.size(); ++i) {
		current += buckets[i];
		if (current > toRetain) {
			cutoff = cutoffs[i];
			break;
		}
	}

	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = 0; j < matrices.size(); ++j) {

			if (i == j) {
				continue;
			}
			auto matrix = dynamic_cast<SparseMatrixType*>(matrices[i][j]);
		//	matrix->filter(cutoff);
		}
	}

	return cutoff;
}

void SparseHelper::filterStatic(Array<SparseMatrixType*>& matrices, float cutoff )
{
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = 0; j < matrices.size(); ++j) {

			if (i == j) {
				continue;
			}
			auto matrix = dynamic_cast<SparseMatrixType*>(matrices[i][j]);
		//	matrix->filter(cutoff);
		}
	}
}

void SparseHelper::getHistogram(SparseMatrixType& matrix, std::map<int,int>& histogram)
{
	for (auto s = matrix.getRowSizes() + 1; s < matrix.getRowSizes() + matrix.getSeq1Length() + 1; ++s) {
		
		bool found;
		for (auto bin = histogram.begin(); bin != histogram.end(); ++bin) {
			if (*s <= bin->first) {
				++histogram[bin->first];
				found = true;
				break;
			} 
		}
		
		if (!found) {
			++histogram[666];
		}	
	}
}

void SparseHelper::getHistogram(Array<SparseMatrixType*>& matrices, std::map<int,int>& histogram)
{
	for (int i = 0; i < matrices.size(); ++i) {
		for (int j = i + 1; j < matrices.size(); ++j) {
			getHistogram(*matrices[i][j], histogram);
		}
	}
}

