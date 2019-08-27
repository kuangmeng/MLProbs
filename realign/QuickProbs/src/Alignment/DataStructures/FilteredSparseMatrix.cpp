#include "FilteredSparseMatrix.h"

#include <cassert>
#include <algorithm>
#include <numeric>

std::vector<float> FilteredSparseMatrix::cutoffs;

void FilteredSparseMatrix::initialise(float initialCutoff, int binsCount)
{
	float step = (1.0f - initialCutoff) / binsCount;

	cutoffs.resize(binsCount);
	
	cutoffs[cutoffs.size() - 1] = initialCutoff;

	for (int i = cutoffs.size() - 1; i > 0; --i) {
		cutoffs[i-1] = cutoffs[i] + step;
	}
}

FilteredSparseMatrix* FilteredSparseMatrix::computeTranspose() const
{
	// create a new sparse matrix
	FilteredSparseMatrix *out = new FilteredSparseMatrix(this->getSeq2Length(), this->getSeq1Length(), this->getNumCells());
	out->fillTransposed(*this);
	return out;
}


/*
int FilteredSparseMatrix::filter(float cutoff)
{
	int remain = 0;
	int maxRow = 0;

	auto input = this->data;
	auto output = this->data;

	// filter output rows
	for (int i = 1; i < height; ++i) {
		int size = rowSizes[i];

		int accepted = 0;
		for (int c = 0; c < size; ++c) {
			if (input->value() >= cutoff) {
				*output = *input;
				++output;
				++accepted;
			}
			++input;
		}

		remain += accepted;
		rowSizes[i] = accepted;

		if (i < height - 1 ) {
			rowIndices[i + 1] = remain;
		}

		if (accepted > maxRow) { maxRow = accepted; } // find longest row
	}

	this->numCells = remain;
	this->rowIndices[0] = remain;
	this->rowSizes[0] = maxRow;

	assert(getMaxRow() < width);
	assert(getNumCells() == std::accumulate(this->rowSizes + 1, this->rowSizes + height, 0));

	return remain;
}


int FilteredSparseMatrix::filterLongRows()
{
	if (getMaxRow() < sparseThreshold) {
		return numCells;
	} else {
		for (int i = 1; i < height; ++i) {
			if (rowSizes[i] > sparseThreshold) {
				++filteredRows;
			}
		}
	}

	// iterate over rows to find ones exceeding threshold
	int maxrow = 0;
	int currentIndex = 0;
	for (int i = 1; i < height; ++i) {
		maxrow = rowSizes[i] > maxrow ? rowSizes[i] : maxrow;

		cell_type *row = data + rowIndices[i]; 
		int toRemove = rowSizes[i] - sparseThreshold;

		// if need copy
		if (currentIndex < rowIndices[i]) {
			std::copy(row, row + rowSizes[i], data + currentIndex);
			rowIndices[i] = currentIndex;
		}

		int fulfil = rowSizes[i];

		if (toRemove > 0) {
			++filteredRows;
			auto minmax = std::minmax_element(row , row + rowSizes[i], [](const cell_type& a, const cell_type& b)->bool {
				return a.second < b.second;
			});
			float hi = (minmax.second)->second; // highest value
			float cur = (minmax.first)->second; // lowest value

			// filter as long as number of retained elements in a row is less than threshold
			do {
				cur = cur * 1.1;
				fulfil = std::count_if(row, row + rowSizes[i], [cur](const cell_type& a)->bool {
					return a.second > cur;
				}); 
			} while (fulfil > sparseThreshold);

			std::remove_if(row, row + rowSizes[i], [cur](const cell_type& a)->bool {
				return a.second > cur;
			});

			rowSizes[i] = fulfil;
		}

		currentIndex += fulfil;
	}

	numCells = currentIndex;
	this->rowIndices[0] = numCells;
	this->rowSizes[0] = maxrow;

	assert(getNumCells() == std::accumulate(this->rowSizes + 1, this->rowSizes + height, 0));
	assert(maxrow < width);

	return numCells;
}
*/


