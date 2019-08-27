#pragma once
#include <vector>

// values written in a backtrack array
#define BACK_ROOT		0
#define BACK_LEFT		1
#define BACK_UP			2
#define BACK_DIAG_2D	3
#define BACK_DIAG_3D	4

class Backtrack
{
public:
	void operator()(
		const std::vector<int>& matrix, 
		int width, 
		int height, 
		int start_x, 
		int start_y,
		std::vector<int>& dotMatrix);


	int calculateMatchCount(
		const std::vector<int>& matrix, 
		char* seq1, 
		char* seq2,
		int width, 
		int height,
		int start_x,
		int start_y);
};