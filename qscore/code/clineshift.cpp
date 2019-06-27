#include "qscore.h"

#if	0
static void DumpRefPair()
	{
	if (g_ptrmsaRef == 0)
		{
		fprintf(stderr, "DumpRefPair: g_ptrmsaRef=NULL\n");
		return;
		}

	const unsigned ColCount = g_ptrmsaRef->GetColCount();
	const unsigned BLOCKLENGTH = 60;
	const unsigned BlockCount = (ColCount + BLOCKLENGTH - 1)/BLOCKLENGTH;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned ColFrom = BlockIndex*BLOCKLENGTH;
		unsigned ColTo = ColFrom + BLOCKLENGTH - 1;
		if (ColTo >= ColCount)
			ColTo = ColCount - 1;
		fprintf(stderr, "%8.8s ", g_ptrmsaRef->GetSeqName(g_RefSeqIndexA));
		for (unsigned ColIndex = ColFrom; ColIndex <= ColTo; ++ColIndex)
			fprintf(stderr, "%c", g_ptrmsaRef->GetChar(g_RefSeqIndexA, ColIndex));
		fprintf(stderr, "\n");

		fprintf(stderr, "%8.8s ", g_ptrmsaRef->GetSeqName(g_RefSeqIndexB));
		for (unsigned ColIndex = ColFrom; ColIndex <= ColTo; ++ColIndex)
			fprintf(stderr, "%c", g_ptrmsaRef->GetChar(g_RefSeqIndexB, ColIndex));
		fprintf(stderr, "\n");
		fprintf(stderr, "\n");
		}
	}
#endif

// Compute the Cline shift score [1] from two pair maps.
// It is the relative simplicity of this code that motivates the use of pair maps.
// [1] Cline, Hughey & Karplus (2002), Bioinformatics 18(2) p.306.
double ClineShift(const int iTestMapA[], const int iRefMapA[], unsigned uLengthA,
  const int iTestMapB[], const int iRefMapB[], unsigned uLengthB, double dEpsilon)
	{
	unsigned uRefPairCount = 0;
	unsigned uTestPairCount = 0;
	double dTotal = 0.0;

	for (unsigned uPosA = 0; uPosA < uLengthA; ++uPosA)
		{
		int iRefPosB = iRefMapA[uPosA];
		if (-1 == iRefPosB)
			continue;

		++uRefPairCount;

		int iTestPosB = iTestMapA[uPosA];
		if (-1 == iTestPosB)
			continue;

		int iShift = iabs(iRefPosB - iTestPosB);
		double dScore = (1 + dEpsilon)/(1 + iShift) - dEpsilon;
		assert(dScore >= -dEpsilon);
		assert(dScore <= 1.0);

		dTotal += dScore;
		}

	for (unsigned uPosB = 0; uPosB < uLengthB; ++uPosB)
		{
		int iTestPosA = iTestMapB[uPosB];
		if (-1 == iTestPosA)
			continue;

		++uTestPairCount;

		int iRefPosA = iRefMapB[uPosB];
		if (-1 == iRefPosA)
			continue;

		int iShift = iabs(iRefPosA - iTestPosA);
		double dScore = (1 + dEpsilon)/(1 + iShift) - dEpsilon;
		assert(dScore >= -dEpsilon);
		assert(dScore <= 1.0);

		dTotal += dScore;
		}

	if (0 == uRefPairCount)
		{
		//DumpRefPair();
		//Quit("ClineShift: No aligned pair in ref alignment");
		return 0.0;
		}

	return dTotal / (double) (uTestPairCount + uRefPairCount);
	}
