// qscore.h

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <errno.h>

#include <algorithm>
#include <vector>
#include <string>

#ifdef _MSC_VER
#include <hash_map>
typedef stdext::hash_map<std::string, unsigned> StrToInt;
#else
#include <ext/hash_map>
#define HASH_MAP	

struct HashStringToUnsigned
	{
	size_t operator()(const std::string &Key)  const
		{
		size_t h = 0;
		size_t Bytes = Key.size();
		for (size_t i = 0; i < Bytes; ++i)
			{
			unsigned char c = (unsigned char) Key[i];
			h = c + (h << 6) + (h << 16) - h;
			}
		return h;
		}
	};

typedef __gnu_cxx::hash_map<std::string, unsigned, HashStringToUnsigned> StrToInt;
#endif

using namespace std;

// Allow different conventions: DEBUG or _DEBUG for debug mode,
// NDEBUG for not debug mode.
#ifdef	_DEBUG
#undef DEBUG
#define	DEBUG	1
#endif

#ifdef	DEBUG
#undef _DEBUG
#define	_DEBUG	1
#endif

#ifdef NDEBUG
#undef	DEBUG
#undef	_DEBUG
#endif

typedef vector<unsigned> IntVec;
typedef vector<bool> BoolVec;

#define	all(t, n)		(t *) allocmem((n)*sizeof(t))
#define	reall(p, t, n)		p = (t *) reallocmem(p, (n)*sizeof(t))
#define zero(p,	t, n)	memset(p, 0, (n)*sizeof(t))
void *allocmem(int bytes);
void freemem(void *p);
void *reallocmem(void *p, int bytes);

static inline bool IsGap(char c)
	{
	return '-' == c || '~' == c || '.' == c || '+' == c || '#' == c;
	}

static inline int iabs(int i)
	{
	return i >= 0 ? i : -i;
	}

class MSA;

const double dInsane = double(0xffffffff);
const unsigned uInsane = 987654321;

unsigned CharToLetter(char c);
char LetterToChar(unsigned Letter);

void ComparePair(const MSA &msaTest, unsigned uTestSeqIndexA,
  unsigned uTestSeqIndexB, const MSA &msaRef, unsigned uRefSeqIndexA,
  unsigned uRefSeqIndexB, double *ptrdSP, double *ptrdPS, double *ptrdCS);

double ComparePairSP(const MSA &msaTest, unsigned uTestSeqIndexA,
  unsigned uTestSeqIndexB, const MSA &msaRef, unsigned uRefSeqIndexA,
  unsigned uRefSeqIndexB);

void ComparePairMap(const int iTestMapA[], const int iTestMapB[],
  const int iRefMapA[], const int iRefMapB[], int iLengthA, int iLengthB,
  double *ptrdSP, double *ptrdPS, double *ptrdCS);

double ComparePairMapSP(const int iTestMapA[], const int iTestMapB[],
  const int iRefMapA[], const int iRefMapB[], int iLengthA, int iLengthB);

double SumPairs(const int iMapRef[], const int iMapTest[], unsigned uLength);

double ClineShift(const int iTestMapA[], const int iRefMapA[], unsigned uLengthA,
  const int iTestMapB[], const int iRefMapB[], unsigned uLengthB, double dEpsilon = 0.2);

void MakePairMaps(const MSA &msaTest, unsigned uTestSeqIndexA, unsigned uTestSeqIndexB,
  const MSA &msaRef, unsigned uRefSeqIndexA, unsigned uRefSeqIndexB, int **ptriTestMapAr,
  int **ptriTestMapBr, int **ptriRefMapAr, int **ptriRefMapBr);

void Quit(const char *Format, ...);
void Warning(const char *Format, ...);

FILE *OpenStdioFile(const char *FileName);
int GetFileSize(FILE *f);

void ParseOptions(int argc, char *argv[]);
bool FlagOpt(const char *Name);
const char *ValueOpt(const char *Name);
const char *RequiredValueOpt(const char *Name);
void Usage();

void CompareMSA(const MSA &msaTest, const MSA &msaRef, double *ptrdSP,
  double *ptrdPS, double *ptrdCS);
double ComputeTC(MSA &msaTest, MSA &msaRef);
void FastQ(const MSA &msaTest, const MSA &msaRef, double &Q, double &TC,
  bool WarnIfNoRefAligned = true);
void ComputeGapScoreMSA(MSA &msaTest, MSA &msaRef, double &GC, double &TC);
void Log(const char *Format, ...);
double PerSeq(const MSA &msaTest, const MSA &msaRef);

void QScore();
void SAB();

#include "msa.h"
#include "seq.h"

extern const char *g_TestFileName;
extern const char *g_RefFileName;

extern const MSA *g_ptrmsaTest;
extern const MSA *g_ptrmsaRef;
extern unsigned g_TestSeqIndexA;
extern unsigned g_TestSeqIndexB;
extern unsigned g_RefSeqIndexA;
extern unsigned g_RefSeqIndexB;
extern bool g_Quiet;
extern bool g_Cline;
extern bool g_SeqDiffWarn;
extern bool g_Verbose;
extern bool g_StripX;
extern bool g_StripZ;
extern bool g_StripB;
