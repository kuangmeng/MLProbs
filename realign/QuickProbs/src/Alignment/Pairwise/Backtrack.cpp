#include "Backtrack.h"


void Backtrack::operator()(
	const std::vector<int>& matrix, 
	int width, 
	int height, 
	int start_x, 
	int start_y,
	std::vector<int>& dotMatrix)
{
	dotMatrix.resize(matrix.size(), 0);

	int x = start_x;
	int y = start_y;
	
	// Second pass is to create alignments and update constrained positions (if any).
	// Matrix indexation is moved with respect to string indexation
	// so strings are accessed after index decrement. Index in alignment is decremented here.
	for (;;) {
		int e = matrix[y * width + x];
		dotMatrix[y * width + x] = 1;

		if (e == BACK_DIAG_2D)		{ --y; --x; } 
		else if (e == BACK_UP)		{ --y; } 
		else if (e == BACK_LEFT)	{ --x; }
		else if (e == BACK_ROOT)	{ break; }
	}
}

int Backtrack::calculateMatchCount(
	const std::vector<int>& matrix, 
	char* seq1, 
	char* seq2, 
	int width, 
	int height, 
	int start_x, 
	int start_y)
{
	int x = start_x;
	int y = start_y;
	int matchCount = 0;

	// Second pass is to create alignments and update constrained positions (if any).
	// Matrix indexation is moved with respect to string indexation
	// so strings are accessed after index decrement. Index in alignment is decremented here.
	for (;;) {
		int e = matrix[y * width + x];
	
		if (e == BACK_DIAG_2D)		{
			--y; --x; 
			if (seq1[y] == seq2[x]) {
				++matchCount;
			}
		} 
		else if (e == BACK_UP)		{ --y; } 
		else if (e == BACK_LEFT)	{ --x; }
		else if (e == BACK_ROOT)	{ break; }
	}

	return matchCount;
}
