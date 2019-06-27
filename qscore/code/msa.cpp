#include "qscore.h"

MSA::MSA()
	{
	m_uSeqCount = 0;
	m_uColCount = 0;
	m_uCacheSeqCount = 0;

	m_szSeqs = 0;
	m_SeqBuffer = 0;
	m_UngapMap = 0;
	m_GapMap = 0;
	m_SeqLengths = 0;
	}

MSA::~MSA()
	{
	Free();
	}

void MSA::Free()
	{
	if (0 != m_UngapMap)
		for (unsigned n = 0; n < m_uSeqCount; ++n)
			delete[] m_UngapMap[n];

	delete[] m_szSeqs;
	delete[] m_SeqBuffer;
	delete[] m_UngapMap;
	m_SeqNames.clear();

	m_uSeqCount = 0;
	m_uColCount = 0;
	m_uCacheSeqCount = 0;

	m_SeqBuffer = 0;
	m_szSeqs = 0;
	}

bool MSA::GetSeqIndex(const char *ptrSeqName, unsigned *ptruSeqIndex) const
	{
	std::map<std::string, unsigned>::const_iterator p =
	  m_SeqNameToIndex.find(ptrSeqName);
	if (p == m_SeqNameToIndex.end())
		return false;
	*ptruSeqIndex = p->second;
	return true;
	}

const char *MSA::GetSeqName(unsigned uSeqIndex) const
	{
#if	_DEBUG
	if (uSeqIndex >= m_uSeqCount)
		Quit("MSA::GetSeqName(%u), count=%u", uSeqIndex, m_uSeqCount);
#endif
	return m_SeqNames[uSeqIndex].c_str();
	}

/***
It is sometimes very convenient to represent a pairwise alignment
as a "pair map", which works as follows.

Let iPos1 be the index into ungapped sequence 1, similarly for iPos2.

Then if a pair of letters (iPos1, iPos2) is aligned:

	iMap1[iPos1] = iPos2 and iMap2[iPos2] = iPos1.

If iPos1 is not in an aligned column, or is aligned to a gap, then
iMap1[iPos1] = -1, and similarly for iMap2. This overloads the meaning
of the integer value, so is questionable software engineering practice;
however it's a simple and convenient solution for small applications.
***/
void MSA::GetPairMap(unsigned uSeqIndex1, unsigned uSeqIndex2, int iMap1[],
  int iMap2[]) const
	{
	assert(uSeqIndex1 < GetSeqCount());
	assert(uSeqIndex2 < GetSeqCount());

	int iPos1 = 0;
	int iPos2 = 0;
	const unsigned uColCount = GetColCount();
	for (unsigned uColIndex = 0; uColIndex < uColCount; ++uColIndex)
		{
		char c1 = GetChar(uSeqIndex1, uColIndex);
		char c2 = GetChar(uSeqIndex2, uColIndex);
		bool bIsGap1 = ::IsGap(c1);
		bool bIsGap2 = ::IsGap(c2);
		if (!bIsGap1 && !bIsGap2)
			{
			if (isupper(c1))
				{
				if (!isupper(c2))
					Quit("Both upper and lower case letters (%c,%c) in ref alignment column %d",
					  c1, c2, uColIndex);
				iMap1[iPos1] = iPos2;
				iMap2[iPos2] = iPos1;
				}
			else
				{
				iMap1[iPos1] = -1;
				iMap2[iPos2] = -1;
				}
			++iPos1;
			++iPos2;
			}
		else if (!bIsGap1 && bIsGap2)
			{
			iMap1[iPos1] = -1;
			++iPos1;
			}
		else if (bIsGap1 && !bIsGap2)
			{
			iMap2[iPos2] = -1;
			++iPos2;
			}
		}

#if	_DEBUG
	{
	int iLength1 = iPos1;
	int iLength2 = iPos2;

	for (int iPos1 = 0; iPos1 < iLength1; ++iPos1)
		{
		int iPos2 = iMap1[iPos1];
		if (-1 == iPos2)
			continue;
		assert(iMap2[iPos2] == iPos1);
		}

	for (int iPos2 = 0; iPos2 < iLength2; ++iPos2)
		{
		int iPos1 = iMap2[iPos2];
		if (-1 == iPos1)
			continue;
		assert(iMap1[iPos1] == iPos2);
		}
	}
#endif
	}

unsigned MSA::GetSeqLength(unsigned uSeqIndex) const
	{
	assert(uSeqIndex < GetSeqCount());
	return m_SeqLengths[uSeqIndex];
	}

bool MSA::IsGap(unsigned uSeqIndex, unsigned uColIndex) const
	{
	return ::IsGap(GetChar(uSeqIndex, uColIndex));
	}

void MSA::SetChar(unsigned uSeqIndex, unsigned uIndex, char c)
	{
#if	_DEBUG
	if (uSeqIndex >= m_uSeqCount || uIndex >= m_uColCount)
		Quit("MSA::GetLetter(%u/%u,%u/%u)", uSeqIndex, m_uSeqCount, uIndex, m_uColCount);
#endif
	m_szSeqs[uSeqIndex][uIndex] = c;
	}

void MSA::MakeUngapMap()
	{
	if (m_UngapMap != 0)
		return;
	m_UngapMap = new unsigned *[m_uSeqCount];
	memset(m_UngapMap, 0, m_uSeqCount*sizeof(unsigned *));
	for (unsigned uSeqIndex = 0; uSeqIndex < m_uSeqCount; ++uSeqIndex)
		MakeUngapMapSeq(uSeqIndex);
	}

void MSA::MakeGapMap()
	{
	if (m_GapMap != 0)
		return;
	m_GapMap = new unsigned *[m_uSeqCount];
	memset(m_GapMap, 0, m_uSeqCount*sizeof(unsigned *));
	for (unsigned uSeqIndex = 0; uSeqIndex < m_uSeqCount; ++uSeqIndex)
		MakeGapMapSeq(uSeqIndex);
	}

void MSA::MakeUngapMapSeq(unsigned uSeqIndex)
	{
	unsigned *ptrMap = new unsigned[m_uColCount];
	memset(ptrMap, 0, m_uColCount*sizeof(unsigned));
	unsigned uUngappedColIndex = 0;
	for (unsigned uGappedColIndex = 0; uGappedColIndex < m_uColCount; ++uGappedColIndex)
		{
		if (IsGap(uSeqIndex, uGappedColIndex))
			ptrMap[uGappedColIndex] = uInsane;
		else
			{
			ptrMap[uGappedColIndex] = uUngappedColIndex;
			++uUngappedColIndex;
			}
		}
	m_UngapMap[uSeqIndex] = ptrMap;
	}

void MSA::MakeGapMapSeq(unsigned uSeqIndex)
	{
	unsigned *ptrMap = new unsigned[m_uColCount];
	memset(ptrMap, 0, m_uColCount*sizeof(unsigned));
	unsigned uUngappedColIndex = 0;
	for (unsigned uGappedColIndex = 0; uGappedColIndex < m_uColCount; ++uGappedColIndex)
		{
		if (!IsGap(uSeqIndex, uGappedColIndex))
			{
			ptrMap[uUngappedColIndex] = uGappedColIndex;
			++uUngappedColIndex;
			}
		}
	m_GapMap[uSeqIndex] = ptrMap;
	}

unsigned MSA::GetUngappedColIndex(unsigned uSeqIndex, unsigned uColIndex)
	{
	assert(uSeqIndex < m_uSeqCount);
	assert(uColIndex < m_uColCount);
	return m_UngapMap[uSeqIndex][uColIndex];
	}

unsigned MSA::GetGappedColIndex(unsigned uSeqIndex, unsigned uUngappedColIndex)
	{
	assert(uSeqIndex < m_uSeqCount);
	assert(uUngappedColIndex < m_uColCount);
	return m_GapMap[uSeqIndex][uUngappedColIndex];
	}
