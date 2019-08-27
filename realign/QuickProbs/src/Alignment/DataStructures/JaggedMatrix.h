#pragma once
#include <vector>
#include <cassert>
#include <cstdio>

#include "Common/mathex.h"

inline 
int jaggedIndex(int x, int y, int height, int jag)
{
	int local_x = x % jag;
	int local_offset = (local_x + y) * jag + local_x;
	int global_offset = (x / jag) * (height + jag - 1) * jag;

	return global_offset + local_offset;
}


template <class T>
class JaggedMatrix
{
public:
	int width;

	int height;

	int jag;

	int blockSize;

	const ::size_t size() const { return storage.size(); }

	T* data() { return storage.data(); }

	const T* data() const { return storage.data(); }

	static const ::size_t size(int width, int height, int jag)
	{
		return (height + jag - 1) * mathex::ceildiv(width, jag) * jag;
	}

	int unjaggedToJagged(int u)
	{
		int x = u % width;
		int y = u / width;

		int offset = (x / jag) * blockSize;
		int local_x = x % jag;
		int local = (local_x + y) * jag + local_x;
		int final = offset + local;

		return final;
	}

	int jaggedToUnjagged(int j)
	{
		int offset = (j / blockSize) * blockSize;

		int local = j % blockSize;
		int local_x = local % jag;

		int y = (local - local_x) / jag - local_x;
		int x =  jag * (j / blockSize) + local_x;

		return y * width + x;
	}

	JaggedMatrix(int width, int height, int jag) 
	: width(width), height(height), jag(jag) 
	{
		int jaggedHeight = height + jag - 1;
		int jaggedWidth = mathex::ceildiv(width, jag) * jag;
		blockSize = jaggedHeight * jag;
		storage.resize(jaggedHeight * jaggedWidth, -1);
	}

	JaggedMatrix(const std::vector<T>& matrix, int width, int height, int jag)
		: width(width), height(height), jag(jag)
	{

		int jaggedHeight = height + jag - 1;
		int jaggedWidth = mathex::ceildiv(width, jag) * jag;
		blockSize = jaggedHeight * jag;

		storage.resize(jaggedHeight * jaggedWidth, -1);

		for (int u = 0; u < width * height; ++u)
		{
			int j = unjaggedToJagged(u);
			int uu = jaggedToUnjagged(j);

			assert(u == uu);

			storage[j] = matrix[u];
		}
	}

	std::vector<T> unpack()
	{
		std::vector<T> out(width * height, 0);

		for (int u = 0; u < width * height; ++u)
		{
			int j = unjaggedToJagged(u);
			out[u] = storage[j];
		}

		return out;
	}

	std::vector<T> storage;
};
