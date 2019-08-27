#pragma once
#include <algorithm>

class BufferSet {
public:
	double* d01()	{ return rawPointer; }
	float* f0()		{ return buffer; }
	float* f1()		{ return buffer + layerSize; }
	float* f2()		{ return buffer + 2 * layerSize; }

	BufferSet(size_t layerSize) : layerSize(layerSize) {
        layerSize += layerSize % 2 + 16; // for alignment reasons
		rawPointer = new double[3 * layerSize / 2];
		buffer = (float*)rawPointer;
	}

	~BufferSet() {
		delete [] rawPointer;
	}

private:
	double* rawPointer;
	float* buffer;
	size_t layerSize;
};
