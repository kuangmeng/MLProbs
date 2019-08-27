#pragma once
#include <CL/cl.h>
#include "Common/mathex.h"

#define jagSize_float	32
#define jagSize_double	16

namespace quickprobs 
{

struct PosteriorSchedule
{
	cl_int taskCount;
};

struct RelaxationSchedule
{
	cl_int numSeqs;
	cl_int sparseWidth;
	cl_int sequenceBegin;
	cl_int sequenceEnd;
};

struct PosteriorTask
{
	cl_uint i;	// first sequence index
	cl_uint j;	// second sequence index

	union {
		cl_uint offset_Aij;	// auxiliary Aij matrix offset
		cl_uint layerSize;	// size of single layer size
	};
	union {
		cl_uint offset_Pij;	// posterior Pij matrix offset
		cl_uint width;			// matrix width
	};
	union {
		cl_uint numCells;	// size of sparse matrix
		cl_uint height;	// matrix height
	};

	cl_float pid;

	
	static ::size_t probabilisticElements(int width)					{ return 3 * mathex::ceilround(width, 32); }
	static ::size_t probabilisticForwardElements(int width)				{ return 3 * mathex::ceilround(width, 32); }
	static ::size_t probabilisticBackwardElements(int width)			{ return 2 * mathex::ceilround(width, 32); }
	static ::size_t partitionElements(int width, bool useDoubles)		{ return (useDoubles ? 8 : 4) * mathex::ceilround(width, 32); }
	static ::size_t partitionForwardElements(int width, bool useDoubles){ return (useDoubles ? 6 : 3) * mathex::ceilround(width, 32); }
	static ::size_t partitionReverseElements(int width, bool useDoubles){ return (useDoubles ? 8 : 4) * mathex::ceilround(width, 32); }
	static ::size_t finalizationElements(int width)						{ return 1 * mathex::ceilround(width, 32); }
	
	static ::size_t localElements(int width, bool useDoubles)			{  return (useDoubles ? 8 : 3) * mathex::ceilround(width, 32); }
	static ::size_t elementsPerSymbol(bool useDoubles)					{  return (useDoubles ? 8 : 3); }
};

struct RelaxationTask
{
	cl_ushort i;				// first sequence index
	cl_ushort j;				// second sequence index
	cl_ushort seed;				// seed for random generator
	cl_ushort acceptedCount;	// number of sequences accepted by selectivity

	static ::size_t localBytes(int maxSparseRow, int stripeCount, int stripeLength);
	static int countMaxSparseLength(int localBytes, int stripeCount, int stripeLength);
};

};