/// <summary>
/// The file facilitates sparse matrices processing
/// </summary>

typedef int index_type;

typedef int column_type;
typedef float value_type;
typedef struct cell_type{
	column_type column;
	value_type value;
};

/*
typedef int relax_column_type;
typedef float relax_value_type;
typedef struct {
	relax_column_type column;
	relax_value_type value;
} relax_cell_type;
typedef cell_type relax_cell_type; 

*/

// basic accessors
global float* sparse_numCells(global float* ptr) {
	return ptr + 1;
}

global index_type* sparse_rowSizes(const global float * const ptr) {
	return (global index_type*)(ptr + 2);
}

global index_type* sparse_rowIndices(const global float * const ptr, const int height) {
	return (global index_type*)(ptr + 2) + height;
}

global struct cell_type* sparse_data(const global float * const ptr, const int height) {
	return (global struct cell_type*)((global index_type*)(ptr + 2) + 2 * height);
}

global struct cell_type* sparse_cell(const global float * const ptr, const int height, const int rowId, const int laneId) {
	return sparse_data(ptr, height) + sparse_rowIndices(ptr, height)[rowId] + laneId;
}

// fixed point accessors
global ushort2* sparse_data_fixed(const global float * const ptr, const int height) {
	return (global ushort2*)((global index_type*)(ptr + 2) + 2 * height);
}

global ushort2* sparse_cell_fixed(const global float * const ptr, const int height, const int rowId, const int laneId) {
	return sparse_data_fixed(ptr, height) + sparse_rowIndices(ptr, height)[rowId] + laneId;
}

// offset accessors
global index_type* sparse_rowSizes_offset(const global float * const ptr) {
	return (global index_type*)(ptr);
}

global index_type* sparse_rowIndices_offset(const global float * const ptr, const int height) {
	return (global index_type*)(ptr) + height;
}

global struct cell_type* sparse_data_offset(const global float * const ptr, const int height) {
	return (global struct cell_type*)((global index_type*)(ptr) + 2 * height);
}

global struct cell_type* sparse_cell_offset(const global float * const ptr, const int height, const int rowId, const int laneId) {
	return sparse_data_offset(ptr, height) + sparse_rowIndices_offset(ptr, height)[rowId] + laneId;
}


/// <summary>
/// </summary>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=posterior></param>
/// <param name=auxiliaryBuffer></param>
/// <param name=sparse></param>
__kernel void SparseMatrix_Generate_long(	
					int seq1Length, 
					int seq2Length,
	global	const	float* posterior,
	local			int* verticalBuffer,
	global			float* sparse)
{
	#define width (seq2Length + 1)
	#define height (seq1Length + 1)
	int j_offset = get_local_id(0);
	
	global index_type* out_rowSizes = sparse_rowSizes(sparse);	
	global index_type* out_rowIndices = sparse_rowIndices(sparse, seq1Length + 1);
	global struct cell_type* out_data =  sparse_data(sparse, seq1Length + 1);

	int numCells = 0;
	int maxRowSize = 0;
	int j = j_offset;
	int lanesCount = _ceildiv(width, get_local_size(0));
	
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
	
		int index = jaggedIndex(j, -j_offset + 1, seq1Length + 1, jagSize_float);			// seek to -j'th row
		global const float *p = posterior + index;

		// calculate number of non-zero elements in each row
		// vertical buffer stores size of rows
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		--iterCount; // omit 0'th column

		for (int iter = 0; iter < iterCount; ++iter)
		{
			int i = iter - j_offset + 1; 
			if (i > 0 && i <= seq1Length)
			{
				if (j == 0)	{ verticalBuffer[i] = 0; }
				else if (j <= seq2Length && *p >= POSTERIOR_CUTOFF) { ++verticalBuffer[i];  }
			}
		
			p += jagSize_float;
			barrier(CLK_GLOBAL_MEM_FENCE);
		}
	}

	// copy verticalBuffer to rowSizes
	for (int q = j_offset; q <= seq1Length; q += get_local_size(0)) {
		out_rowSizes[q] = verticalBuffer[q]; 
	}
	barrier(CLK_GLOBAL_MEM_FENCE);

	// calculate starting indices in sparse matrix 
	if (j_offset == 0) {
		verticalBuffer[0] = 0;
		
		for (int i = 0; i <= seq1Length; ++i) {
			int size = verticalBuffer[i];
			maxRowSize = max(maxRowSize, size);
			verticalBuffer[i] = numCells;
			numCells += size;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// copy verticalBuffer to rowIndices
	for (int q = j_offset; q <= seq1Length; q += get_local_size(0)) {
		out_rowIndices[q] = verticalBuffer[q]; 
	}
	barrier(CLK_GLOBAL_MEM_FENCE);

	j = j_offset;
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
	
		int index = jaggedIndex(j, -j_offset + 1, seq1Length + 1, jagSize_float);			// seek to -j'th row
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		--iterCount; // omit 0'th column

		// fill in data matrix
		global const float* p = posterior + index;
		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset + 1; 
			if (i > 0 && i <= seq1Length && j > 0 && j <= seq2Length) {
				float v = *p;
				if (v >= POSTERIOR_CUTOFF) {
					out_data[verticalBuffer[i]].column = j;
					out_data[verticalBuffer[i]].value = v;
					++verticalBuffer[i];
				} 
			}
			p += jagSize_float;
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		}
	}
	
	if (j_offset == 0) {
		out_rowIndices[0] = numCells; // store number of non-zero cells here
		out_rowSizes[0] = maxRowSize;	// store maximum row size

		sparse[0] = posterior[0];	// store distances
		sparse[1] = numCells;		// store number of non-zero elements
	}

	barrier(CLK_GLOBAL_MEM_FENCE);
}

/// <summary>
/// </summary>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=posterior></param>
/// <param name=auxiliaryBuffer></param>
/// <param name=sparse></param>
__kernel void SparseMatrix_Generate_verylong(	
					int seq1Length, 
					int seq2Length,
	global	const	float* posterior,
	local			int* verticalBuffer,
	global			float* sparse)
{
	#define width (seq2Length + 1)
	#define height (seq1Length + 1)
	int j_offset = get_local_id(0);
	
	global index_type* out_rowSizes = sparse_rowSizes(sparse);	
	global index_type* out_rowIndices = sparse_rowIndices(sparse, seq1Length + 1);
	global struct cell_type* out_data =  sparse_data(sparse, seq1Length + 1);

	int numCells = 0;
	int maxRowSize = 0;
	int j = j_offset;
	int lanesCount = _ceildiv(width, get_local_size(0));
	
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
	
		int index = jaggedIndex(j, -j_offset + 1, seq1Length + 1, jagSize_float);			// seek to -j'th row
		global const float *p = posterior + index;

		// calculate number of non-zero elements in each row
		// vertical buffer stores size of rows
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		--iterCount; // omit 0'th column

		for (int iter = 0; iter < iterCount; ++iter)
		{
			int i = iter - j_offset + 1; 
			if (i > 0 && i <= seq1Length)
			{
				if (j == 0)	{ out_rowSizes[i] = 0; }
				else if (j <= seq2Length && *p >= POSTERIOR_CUTOFF) { ++out_rowSizes[i];  }
			}
		
			p += jagSize_float;
			barrier(CLK_GLOBAL_MEM_FENCE);
		}
	}

	// calculate starting indices in sparse matrix 
	if (j_offset == 0) {
		out_rowIndices[0] = 0;
		out_rowSizes[0] = 0;
		
		for (int i = 0; i <= seq1Length; ++i) {
			int size = out_rowSizes[i];
			maxRowSize = max(maxRowSize, size);
			out_rowIndices[i] = numCells;
			numCells += size;
		}
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	j = j_offset;
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
	
		int index = jaggedIndex(j, -j_offset + 1, seq1Length + 1, jagSize_float);			// seek to -j'th row
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		--iterCount; // omit 0'th column

		// fill in data matrix
		global const float* p = posterior + index;
		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset + 1; 
			if (i > 0 && i <= seq1Length && j > 0 && j <= seq2Length) {
				float v = *p;
				if (v >= POSTERIOR_CUTOFF) {
					out_data[out_rowIndices[i]].column = j;
					out_data[out_rowIndices[i]].value = v;
					++out_rowIndices[i];
				} 
			}
			p += jagSize_float;
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		}
	}
	
	// restore offsets
	for (int i = j_offset; i <= seq1Length; i += get_local_size(0)) {
		out_rowIndices[i] -= out_rowSizes[i];
	}

	barrier(CLK_GLOBAL_MEM_FENCE);

	if (j_offset == 0) {
		out_rowIndices[0] = numCells; // store number of non-zero cells here
		out_rowSizes[0] = maxRowSize;	// store maximum row size

		sparse[0] = posterior[0];	// store distances
		sparse[1] = numCells;		// store number of non-zero elements
	}

	barrier(CLK_GLOBAL_MEM_FENCE);
}

/// <summary>
/// </summary>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=posterior></param>
/// <param name=auxiliaryBuffer></param>
/// <param name=sparse></param>
__kernel void SparseMatrix_GenerateTranspose_long(	
					int seq1Length, 
					int seq2Length,
	global	const	float* posterior,
	local			int* horizontalBuffer,
	global			float* sparse)
{
	int j_offset = get_local_id(0);
	
	global index_type* out_rowSizes = sparse_rowSizes(sparse);	
	global index_type* out_rowIndices = sparse_rowIndices(sparse, seq2Length + 1);
	global struct cell_type* out_data =  sparse_data(sparse, seq2Length + 1);

	int numCells = 0;
	int maxRowSize = 0;
	
	int j = j_offset;
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
	
		int rowSize = 0;
		int index = jaggedIndex(j, -j_offset + 1, seq1Length + 1, jagSize_float);		// omit 0'th row
		global const float *p = posterior + index;

		// calculate number of non-zero elements in each column
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		--iterCount; // omit 0'th column

		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset + 1; 
			if (i > 0 && i <= seq1Length) {
				if (j > 0 && j <= seq2Length && *p >= POSTERIOR_CUTOFF) { ++rowSize; } 
			}
		
			p += jagSize_float;
			barrier(CLK_GLOBAL_MEM_FENCE);
		}
	
		// copy register value to rowSizes and horizontalBuffer
		if (j <= seq2Length) { out_rowSizes[j] = horizontalBuffer[j_offset] = rowSize; }
		barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

		// calculate starting indices in sparse matrix
		if (j_offset == 0) {
			int end = min((int)get_local_size(0), seq2Length + 1 - lane * (int)get_local_size(0));
			for (int x = 0; x < end; ++x) {
				int size = horizontalBuffer[x];
				maxRowSize = max(maxRowSize, size);
				horizontalBuffer[x] = numCells;
				numCells += size;
			}
		}

		barrier(CLK_LOCAL_MEM_FENCE);
	
		// copy horizontal to rowIndices
		if (j <= seq2Length) { out_rowIndices[j] = horizontalBuffer[j_offset]; }
		barrier(CLK_GLOBAL_MEM_FENCE);

		// fill in data matrix
		p = posterior + index;
		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset + 1; 
			if (i > 0 && i <= seq1Length && j > 0 && j <= seq2Length) {
				float v = *p;
				if (v >= POSTERIOR_CUTOFF) {
					out_data[horizontalBuffer[j_offset]].column = i;
					out_data[horizontalBuffer[j_offset]].value = v;
					++horizontalBuffer[j_offset];
				} 
			}
			p += jagSize_float;
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		}
	}

	if (j_offset == 0) {
		out_rowIndices[0] = numCells; // store number of non-zero cells here
		out_rowSizes[0] = maxRowSize;	// store maximum row size

		sparse[0] = posterior[0];	// store distances
		sparse[1] = numCells;		// store number of non-zero elements
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
}

/*
/// <summary>
/// </summary>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=data></param>
/// <param name=rowSize></param>
/// <param name=posterior></param>
/// <returns></returns>
__kernel void SparseMatrix_ToDense(	
					int seq1Length, 
					int seq2Length, 
	global	const	struct cell_type* data,
	global	const	int* rowSizes,
	global			float* posterior)
{	
	int j = get_local_id(0);
	int index = jaggedIndex(j, -j, seq1Length + 1, jagSize_float);			// seek to -j'th row
	
	int iterationsCount  = seq1Length + seq2Length + 1;
	for (int iter = 0; iter < iterationsCount; ++iter) {
		int i = iter - j;
		
		if (i >= 0 && i <= seq1Length && j <= seq2Length) {
			posterior[index] = 0;
		}

		index += jagSize_float;
	}

	barrier(CLK_GLOBAL_MEM_FENCE);
	
	global struct cell_type* rowPtr = data;
	for (int i = 1; i <= seq1Length; i++) {
		if (j < rowSizes[i]) {
			int index = jaggedIndex(rowPtr[j].column, i, seq1Length + 1, jagSize_float);
			posterior[index] = rowPtr[j].value;
		}
		
		rowPtr += rowSizes[i];
	}
}
*/