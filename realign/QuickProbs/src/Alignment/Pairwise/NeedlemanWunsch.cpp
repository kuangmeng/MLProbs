#include <vector>
#include <limits>

#include "Backtrack.h"
#include "NeedlemanWunsch.h"
#include "Common/mathex.h"

/// <summary>
/// See declaration for all the details.
/// </summary>
int NeedlemanWunsch::operator()(						
	const char *seq1, 
	const char *seq2,
	int seq1Length,
	int seq2Length)
{
	int w = seq2Length + 1;
	int h = seq1Length + 1;
	
	int *auxBuffer = new int[2 * w];
	int *F = auxBuffer;
	int *H = auxBuffer + w;


	// reset first row and first column of all arrays,
	// note row-major addressing
	F[0] = 0;
	H[0] = 0;
	
	for (int x = 1; x < w; x++)
		F[x] = H[x] = (x - 1) * ge + gi;	// H[0] = 0, H[1] = gi, H[2] = gi + ge

	// main dynamic programming recursion
	for (int y = 1; y < h; y++) {
		int diagH = H[0];			// H[0] from previous iteration is diagonal
		H[0] = (y - 1) * ge + gi;	// recalculate H[0]
		int e = (y - 1) * ge + gi;

		for (int x = 1; x < w; x++) {
			// before recalculation F[x], H[x] represent elements from upper row
			e = mathex::max(e + ge, H[x - 1] + gi);
			F[x] = mathex::max(F[x] + ge, H[x] + gi);
			int diag2d = diagH + substitution[seq1[y - 1] * symbolsCount + seq2[x - 1]];	

			diagH = H[x]; // current from will be diagonal after H[x] recalculation and x increment
			H[x] = mathex::max(e, F[x], diag2d);
		}
	}

	int score = H[w-1];
	delete [] auxBuffer;

	return score;
}

/// <summary>
/// See declaration for all the details.
/// </summary>
int NeedlemanWunsch::operator()(
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

	// reset first row and first column of all arrays (including backtrack array) 
	E[0] = F[0] = H[0] = 0;
	backtrack[0] = BACK_ROOT;

	for (int x = 1; x < w; x++) {
		E[x] = F[x] = H[x] = (x - 1) * ge + gi;
		backtrack[x] = BACK_LEFT;
	}

	for (int y = 1; y < h; y++) {
		E[y * w] = F[y * w] = H[y * w] = (y - 1) * ge + gi;
		backtrack[y * w] = BACK_UP;
	}

	// main dynamic programming recursion
	for (int y = 1; y < h; y++)	{
		int idx = y * w + 1;
		int idxUp = (y - 1) * w + 1;

		for (int x = 1; x < w; x++) {
			int w = substitution[seq1[y - 1] * symbolsCount +seq2[x - 1]];

			E[idx] = mathex::max(E[idx - 1] + ge, H[idx - 1] + gi);	// left variant
			F[idx] = mathex::max(F[idxUp] + ge, H[idxUp] + gi);		// up variant

			H[idx] = mathex::max3idx(E[idx], F[idx], H[idxUp - 1] + w, backtrack[idx]);
			backtrack[idx]++;

			idx++;
			idxUp++;
		}
	}

	// backtrack in order to obtain the best alignment

	int score = H[h * w - 1];
	delete [] auxBuffer;
	
	return score;
}