#include "qscore.h"

#define UINT_MAX   4294967295U

unsigned g_SeqDiffCount;

// O(NL) computation of PREFAB Q score and Balibase TC score.
// Algorithm based on an idea due to Chuong (Tom) Do.
// Each position in the reference alignment is annotated with
// the column number C in the test alignment where the same
// letter is found. A pair of identical Cs in the same reference
// column indicates a correctly aligned pair of letters.
void FastQ(const MSA &msaTest, const MSA &msaRef, double &Q, double &TC,
  bool WarnIfNoRefAligned)
	{
	g_SeqDiffCount = 0;

	unsigned CorrectPairCount = 0;
	unsigned RefAlignedPairCount = 0;

	const unsigned RefSeqCount = msaRef.GetSeqCount();
	const unsigned TestSeqCount = msaTest.GetSeqCount();

	const unsigned RefColCount = msaRef.GetColCount();
	const unsigned TestColCount = msaTest.GetColCount();

	StrToInt RefSeqNameToIndex;
	IntVec RefToTestSeqIndex(RefSeqCount);

	for (unsigned RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		const string SeqName = msaRef.GetSeqName(RefSeqIndex);
		RefToTestSeqIndex[RefSeqIndex] = UINT_MAX;
		RefSeqNameToIndex[SeqName] = RefSeqIndex;
		}

	for (unsigned TestSeqIndex = 0; TestSeqIndex < TestSeqCount; ++TestSeqIndex)
		{
		const string SeqName = msaTest.GetSeqName(TestSeqIndex);
		StrToInt::const_iterator p =
		  RefSeqNameToIndex.find(SeqName);
		if (p != RefSeqNameToIndex.end())
			{
			unsigned RefSeqIndex = p->second;
			if (RefSeqIndex == UINT_MAX)
				Quit("UINT_MAX");
			RefToTestSeqIndex[RefSeqIndex] = TestSeqIndex;
			}
		}

#ifdef DEBUG
	if (!FlagOpt("ignoremissingseqs"))
	for (unsigned RefSeqIndex = 0; RefSeqIndex < RefSeqCount; ++RefSeqIndex)
		{
		unsigned TestSeqIndex = RefToTestSeqIndex[RefSeqIndex];
		const char *RefSeqName = msaRef.GetSeqName(RefSeqIndex);
		if (TestSeqIndex == UINT_MAX)
			{
			fprintf(stderr, "\n");
			fprintf(stderr, "RefSeqIndex  RefSeqName\n");
			fprintf(stderr, "===========  ==========\n");
			for (StrToInt::const_iterator p =
			  RefSeqNameToIndex.begin(); p != RefSeqNameToIndex.end(); ++p)
				fprintf(stderr, "%11u  %s\n", p->second, (p->first).c_str());
			fprintf(stderr, "\n");
			Quit("Ref seq %u=%.16s not found in test alignment", RefSeqIndex, RefSeqName);
			}
		else
			{
			const char *TestSeqName = msaTest.GetSeqName(TestSeqIndex);
			assert(!strcmp(RefSeqName, TestSeqName));
			}
		}
#endif

// TestColIndex[i] is the one-based (not zero-based!) test column index
// of the letter found in the current column of the reference alignment
// (or the most recent letter if the reference column is gapped, or zero
// if no letter has yet been found). Here, seq index i is for msaRef.
	IntVec TestColIndex(TestSeqCount, 0);

// TestColIndexCount[i] is the number of times that a letter from test
// column i (one-based!) appears in the current reference column.
	IntVec TestColIndexCount(TestColCount+1, 0);

// TestColIndexes[i] is the column index in the test alignment of
// the i'th non-gapped position in the current reference column.
	IntVec TestColIndexes;

	unsigned RefAlignedColCount = 0;
	unsigned CorrectColCount = 0;
	if (g_Verbose)
		{
		Log("RefCol  RefAln  NonGapped  TestAll  CorrCols  Ref\n");
		Log("------  ------  ---------  -------  --------  ---\n");
		//        6          9        7         8
		}

	for (unsigned RefColIndex = 0; RefColIndex < RefColCount; RefColIndex++)
		{
		TestColIndexes.clear();
		TestColIndexes.reserve(RefSeqCount);

	// NonGappedCount is the number of non-gapped positions in the current
	// reference column.
		unsigned NonGappedCount = 0;
		unsigned FirstTestColIndex = UINT_MAX;
		bool RefColIsAligned = false;
		bool TestColAllCorrect = true;
		bool TestAllAligned = true;
		for (unsigned RefSeqIndex = 0; RefSeqIndex < RefSeqCount; RefSeqIndex++)
			{
			unsigned TestSeqIndex = RefToTestSeqIndex[RefSeqIndex];
			if (TestSeqIndex == UINT_MAX)
				{
				if (FlagOpt("ignoremissingseqs"))
					continue;
				Quit("Test seq %.16s missing", msaRef.GetSeqName(RefSeqIndex));
				}

			char cRef = msaRef.GetChar(RefSeqIndex, RefColIndex);
			if (!IsGap(cRef))
				{
				char cTest = 0;
				unsigned Col = TestColIndex[TestSeqIndex];
				do
					cTest = msaTest.GetChar(TestSeqIndex, Col++);
				while (IsGap(cTest));
				if (toupper(cRef) != toupper(cTest))
					{
					if (g_SeqDiffWarn)
						{
						++g_SeqDiffCount;
#if 0
						Warning("Test seq %u (%s) differs from ref seq %u (%s), ref col %u=%c, test=%c",
						  TestSeqIndex,
						  msaTest.GetSeqName(TestSeqIndex),
						  RefSeqIndex,
						  msaRef.GetSeqName(RefSeqIndex),
						  RefColIndex,
						  cRef,
						  cTest);
#endif
						}
					else
						Quit("Test seq %u (%s) differs from ref seq %u (%s), ref col %u=%c, test=%c",
						  TestSeqIndex,
						  msaTest.GetSeqName(TestSeqIndex),
						  RefSeqIndex,
						  msaRef.GetSeqName(RefSeqIndex),
						  RefColIndex,
						  cRef,
						  cTest);
					}
				if (isalpha(cRef) && isupper(cRef))
					{
					RefColIsAligned = true;
					++NonGappedCount;
					if (isupper(cTest))
						{
						TestColIndexes.push_back(Col);
						++(TestColIndexCount[Col]);
						if (FirstTestColIndex == UINT_MAX)
							FirstTestColIndex = Col;
						else
							{
							if (FirstTestColIndex != Col)
								TestColAllCorrect = false;
							}
						}
					else
						TestAllAligned = false;
					}
				else
					{
					if (RefColIsAligned)
						{
						fprintf(stderr, "\n");
						fprintf(stderr, "Ref col: ");
						for (unsigned RefSeqIndex = 0; RefSeqIndex < RefSeqCount; RefSeqIndex++)
							fprintf(stderr, "%c", msaRef.GetChar(RefSeqIndex, RefColIndex));
						fprintf(stderr, "\n");
						Quit("Ref col %u has both upper- and lower-case letters",
						  RefColIndex);
						}
					}
				TestColIndex[TestSeqIndex] = Col;
				}
			}

		if (RefColIsAligned && NonGappedCount > 1)
			{
			++RefAlignedColCount;
			if (TestColAllCorrect && TestAllAligned)
				++CorrectColCount;
			}

		unsigned ColPairCount = 0;
		for (IntVec::const_iterator p = TestColIndexes.begin(); p != TestColIndexes.end(); ++p)
			{
			unsigned Col = *p;
			unsigned Count = TestColIndexCount[Col];
			if (Count > 0)
				ColPairCount += Count*(Count - 1)/2;
			TestColIndexCount[Col] = 0;
			}

		CorrectPairCount += ColPairCount;
		RefAlignedPairCount += NonGappedCount*(NonGappedCount - 1)/2;

		if (g_Verbose)
			{
			Log("%6u  %6c  %9u  %7c  %8u  ",
			  RefColIndex,
			  RefColIsAligned ? 'T' : 'F',
			  NonGappedCount,
			  TestColAllCorrect ? 'T' : 'F',
			  CorrectColCount);
			for (unsigned RefSeqIndex = 0; RefSeqIndex < RefSeqCount; RefSeqIndex++)
				{
				unsigned TestSeqIndex = RefToTestSeqIndex[RefSeqIndex];
				if (TestSeqIndex == UINT_MAX)
					continue;

				char cRef = msaRef.GetChar(RefSeqIndex, RefColIndex);
				Log("%c", cRef);
				}
			Log("\n");
			}
		}

	if (RefAlignedPairCount == 0)
		Q = 0;
	else
		Q = (double) CorrectPairCount / (double) RefAlignedPairCount;

	if (RefAlignedColCount == 0)
		{
		if (WarnIfNoRefAligned)
			fprintf(stderr,
			  "Warning: reference alignment %s has no aligned (upper-case) columns\n",
			  g_RefFileName);
		TC = 0;
		}
	else
		TC = (double) CorrectColCount / (double) RefAlignedColCount;

	if (g_Verbose)
		{
		Log("        ------                      --------\n");
		Log("%6.6s  %6u  %9.9s  %7.7s  %8u\n",
		  "", RefAlignedColCount, "", "", CorrectColCount);
		Log("\n");
		Log("CorrectPairCount     %u\n", CorrectPairCount);
		Log("RefAlignedPairCount  %u\n", RefAlignedPairCount);
		Log("CorrectColCount      %u\n", CorrectColCount);
		Log("RefAlignedColCount   %u\n", RefAlignedColCount);
		Log("Q                    %.4f\n", Q);
		Log("TC                   %.4f\n", TC);
		}
	if (g_SeqDiffCount > 0)
		Warning("%u seq diffs ignored", g_SeqDiffCount);
	}
