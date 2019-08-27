#pragma once

// some stuff related to global memory pattern
#define jagSize_float	32
#define jagSize_double	16

#define JAG_FLOAT		32
#define JAG_DOUBLE		16


int jaggedSize(int width, int height, int jag)
{
	return (height + jag - 1) * _ceildiv(width, jag) * jag;
}

int jaggedIndex(int x, int y, int height, int jag)
{
	int local_x = x % jag;
	int local_offset = (local_x + y) * jag + local_x;
	int global_offset = (x / jag) * (height + jag - 1) * jag;
	
	return global_offset + local_offset;
}

#define jaggedIndex_float(x,y,height) jaggedIndex(x,y,height, jagSize_float)

/*
int jaggedIndex_float(int x, int y, int height)
{
	int local_x = x & 31;
	int local_offset = (local_x + y) * 32 + local_x;
	int global_offset = (x >> 5) * (height + 31) << 5;
	
	return global_offset + local_offset;
}
*/

struct PosteriorSchedule
{
	int taskCount;
};


struct RelaxationSchedule
{
	int numSeqs;
	int sparseWidth;
	int sequenceBegin;
	int sequenceEnd;
};

struct PosteriorTask
{
	unsigned int i;	// first sequence index
	unsigned int j;	// second sequence index
	unsigned int offset_Aij;	// auxiliary buffer offset
	unsigned int offset_Pij;	// posterior Pij matrix offset
	unsigned int numCells;	// size of sparse matrix
	float pid;
};

struct RelaxationTask
{
	unsigned short i;	// first sequence index
	unsigned short j;	// second sequence index
	unsigned short seed;			// random number generator seed
	unsigned short acceptedCount;	// number of sequences accepted by selectivity
};

