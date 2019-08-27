#include "AuxiliaryStructures.h"
#include "DataStructures/SparseMatrixType.h"

using namespace quickprobs;

::size_t RelaxationTask::localBytes(int rowLength, int stripeCount, int stripeLength) 
{
	::size_t bytes = 
		stripeCount * rowLength * sizeof(OpenCLSparseMatrix::cell_type) + // sparse elements of Sxy rows
		stripeCount * stripeLength * sizeof(OpenCLSparseMatrix::cell_type); // sparse elements of Sxz rows

	return bytes;
}

int RelaxationTask::countMaxSparseLength(int localBytes, int stripeCount, int stripeLength)
{
	return localBytes / (sizeof(OpenCLSparseMatrix::cell_type) * stripeCount) - stripeLength - 1;
}