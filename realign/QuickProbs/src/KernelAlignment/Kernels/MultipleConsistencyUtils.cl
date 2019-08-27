

__kernel __attribute__((reqd_work_group_size(STRIPE_COUNT * STRIPE_LENGTH, 1, 1)))
void MultipleConsistency_Postprocess(
global		const	unsigned int* lengths,
global		const	unsigned int* sparseOffsets,
global				struct RelaxationTask* tasks,
global				float* Sxy,								// output sparse matrix to be filtered
			const	float posteriorCutoff,						
					struct RelaxationSchedule schedule,		
local				struct cell_type* Sxy_row)
{
	#define localId get_local_id(0)
	#define taskId get_group_id(0)
	#define stripeId (get_local_id(0) / STRIPE_LENGTH)
	
	int stripeOffset = get_local_id(0) % STRIPE_LENGTH; 
						
	local index_type Sxy_sizes[STRIPE_COUNT];
	local int x, y;
	local int length_x_inc;

	if (localId == 0) {
		x = tasks[taskId].i; 
		y = tasks[taskId].j;
		length_x_inc = lengths[x] + 1;
	} 

	barrier(CLK_LOCAL_MEM_FENCE);

	// copy indices to auxiliary array
	Sxy += sparseOffsets[x * schedule.numSeqs + y];	// output sparse matrix
	
	float sumW = *Sxy;

	int end_ix = _ceilround(length_x_inc - 1, STRIPE_COUNT) + 1;

	// filter rows and update sizes	
	for (int ix = 1 + stripeId; ix < end_ix; ix += STRIPE_COUNT) {	
		
		if (stripeOffset == 0) { 
			Sxy_sizes[stripeId] = ix < length_x_inc ? sparse_rowSizes(Sxy)[ix] : 0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		int Sxy_size = Sxy_sizes[stripeId];

		// store locally Sxy rows being processed and weight them
		for (int y_lane = stripeOffset; y_lane < Sxy_size; y_lane += STRIPE_LENGTH) {
			sparse_cell(Sxy, length_x_inc, ix, y_lane)->value /= sumW;
		}
		
		barrier(CLK_GLOBAL_MEM_FENCE);
	}

	#undef localId
	#undef taskId
	#undef stripeId
}


__kernel __attribute__((reqd_work_group_size(STRIPE_COUNT * STRIPE_LENGTH, 1, 1)))
void MultipleConsistency_PostprocessFiltered(
global		const	unsigned int* lengths,
global		const	unsigned int* sparseOffsets,
global				struct RelaxationTask* tasks,
global				float* Sxy,								// output sparse matrix to be filtered
			const	float posteriorCutoff,						
					struct RelaxationSchedule schedule,		
local				struct cell_type* Sxy_row)
{
	#define localId get_local_id(0)
	#define taskId get_group_id(0)
	#define stripeId (get_local_id(0) / STRIPE_LENGTH)
	
	int stripeOffset = get_local_id(0) % STRIPE_LENGTH; 
						
	local index_type Sxy_sizes[STRIPE_COUNT];
	local int x, y;
	local int length_x_inc;
	local int totalCells;
	local int stripeRowIndex;
	int maxRow = 0;
	

	if (localId == 0) {
		x = tasks[taskId].i; 
		y = tasks[taskId].j;
		length_x_inc = lengths[x] + 1;
		totalCells = 0;
		stripeRowIndex = 0;
	} 

	barrier(CLK_LOCAL_MEM_FENCE);

	// copy indices to auxiliary array
	Sxy += sparseOffsets[x * schedule.numSeqs + y];	// output sparse matrix
	float sumW = *Sxy;

	Sxy_row += stripeId * schedule.sparseWidth;

	int end_ix = _ceilround(length_x_inc - 1, STRIPE_COUNT) + 1;

	// filter rows and update sizes	
	for (int ix = 1 + stripeId; ix < end_ix; ix += STRIPE_COUNT) {	
		
		if (stripeOffset == 0) { 
			Sxy_sizes[stripeId] = ix < length_x_inc ? sparse_rowSizes(Sxy)[ix] : 0;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		int Sxy_size = Sxy_sizes[stripeId];

		// store locally Sxy rows being processed and weight them
		for (int y_lane = stripeOffset; y_lane < Sxy_size; y_lane += STRIPE_LENGTH) {
			Sxy_row[y_lane] = sparse_data(Sxy, length_x_inc)[
				(stripeId == 0 ? stripeRowIndex : sparse_rowIndices(Sxy, length_x_inc)[ix]) + y_lane]; 
			// use source index to get row as destination indices are altered
			Sxy_row[y_lane].value = Sxy_row[y_lane].value / sumW; 
		}
		
		barrier(CLK_LOCAL_MEM_FENCE);
		
		// filter output rows
		if (stripeOffset == 0 && ix < length_x_inc) {
			int accepted = 0;
			for (int c = 0; c < Sxy_size; ++c) {
				if (Sxy_row[c].value >= posteriorCutoff) {
					Sxy_row[accepted++] = Sxy_row[c];
				}
			}
			
			sparse_rowSizes(Sxy)[ix] = Sxy_sizes[stripeId] = accepted;
		}
		barrier(CLK_LOCAL_MEM_FENCE);
		Sxy_size = Sxy_sizes[stripeId];
		
		// update offsets of Sxy rows being processed
		if (localId == 0) {
			int idx = sparse_rowIndices(Sxy, length_x_inc)[ix];
			stripeRowIndex = (ix + STRIPE_COUNT < length_x_inc) ? sparse_rowIndices(Sxy, length_x_inc)[ix + STRIPE_COUNT] : 0;

			for (int j = 0, row = ix; j < STRIPE_COUNT && row < length_x_inc; ++j, ++row) {
				idx += Sxy_sizes[j];
				maxRow = max(maxRow, Sxy_sizes[j]);
			
				if (row < length_x_inc - 1) {
					sparse_rowIndices(Sxy, length_x_inc)[row + 1] = idx; // update index of next row 
				}
			}
			totalCells = idx;
		}

		barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

		// store output Sxy rows
		for (int y_lane = stripeOffset; y_lane < Sxy_size; y_lane += STRIPE_LENGTH) {
			*sparse_cell(Sxy, length_x_inc, ix, y_lane) = Sxy_row[y_lane]; 
		}

		barrier(CLK_GLOBAL_MEM_FENCE);
	}

	// store new number of cells
	if (localId == 0) {
		*sparse_numCells(Sxy) = totalCells;
		sparse_rowSizes(Sxy)[0] = maxRow; 
		sparse_rowIndices(Sxy, length_x_inc)[0] = totalCells;
	}
	
	#undef localId
	#undef taskId
	#undef stripeId
}

