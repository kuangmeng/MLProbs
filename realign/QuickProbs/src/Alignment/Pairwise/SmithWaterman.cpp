#include <cstdlib>
#include <cstring>
#include <cassert>
#include <limits>

#include "Common/mathex.h"
#include "SmithWaterman.h"

/// <summary>
/// See declaration for all the details.
/// </summary>
int SmithWaterman::operator()(						
	const char *seq1, 
	const char *seq2,
	int seq1Length,
	int seq2Length)
{
	int w = seq2Length + 1;
	int h = seq1Length + 1;

	int *auxBuffer = new int[2 * w];
	int *F = auxBuffer + 0 * w;
	int *H = auxBuffer + 1 * w;

	// reset first row and first column of all arrays,
	// note row-major addressing 
	for (int x = 0; x < w; x++) {
		F[x] = H[x] = 0;
	}

	// main dynamic programming recursion
	int maxH = std::numeric_limits<int>::min();
	
	for (int y = 1; y < h; y++) {
		int e = 0;
		int diagH = 0;
		
		// H[0] is always equal to 0
		for (int x = 1; x < w; x++) {
			// before recalculation F[x], H[x] represent elements from upper row
			e = mathex::max(e + ge, H[x - 1] + gi);
			F[x] = mathex::max(F[x] + ge, H[x] + gi);
			int diag2d = diagH + substitution[seq1[y - 1] * symbolsCount + seq2[x - 1]];	
			
			diagH = H[x]; // current from will be diagonal after H[x] recalculation and x increment
			H[x] = mathex::max(0, e, F[x], diag2d);
			maxH = mathex::max(H[x], maxH);
		}
	}

	delete [] auxBuffer;

	return maxH;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
int SmithWaterman::operator()(
	const char *seq1, 
	const char *seq2,
	int seq1Length,
	int seq2Length,
	std::vector<int>& scores,
	std::vector<int>& backtrack)
{
	int w = seq2Length + 1;
	int h = seq1Length + 1;
	
	scores.resize(w * h);
	backtrack.resize(w * h);

	int *auxBuffer = new int[h * w * 2]; // allocate memory at once
	int *E = auxBuffer + 0 * h * w;
	int *F = auxBuffer + 1 * h * w;
	int *H = scores.data();
	
	start_x = 0;
	start_y = 0;

	// reset first row and first column of all arrays (including backtrack array)
	// note row-major addressing 
	for (int x = 0; x < w; x++) {
		E[x] = F[x] = H[x] = backtrack[x] = 0;
	}

	for (int y = 1; y < h; y++) {
		E[y * w] = F[y * w] = H[y * w] = backtrack[y * w] = 0;
	}

	// main dynamic programming recursion
	int maxH = std::numeric_limits<int>::min();
	
	for (int y = 1; y < h; y++) {
		int idx = y * w + 1;
		int idxUp = (y - 1) * w + 1;

		for (int x = 1; x < w; x++) {
			int w = substitution[seq1[y - 1] * symbolsCount + seq2[x - 1]];
			
			E[idx] = mathex::max(E[idx - 1] + ge, H[idx - 1] + gi);	// left variant
			F[idx] = mathex::max(F[idxUp] + ge, H[idxUp] + gi);		// up variant

			H[idx] = mathex::max4idx(0, E[idx], F[idx],	H[idxUp - 1] + w, backtrack[idx]);
			
			if (H[idx] > maxH) {
				maxH = H[idx];
				start_x = x;
				start_y = y;
			}

			idx++;
			idxUp++;
		}
	}

	delete [] auxBuffer;
	
	return maxH;
}