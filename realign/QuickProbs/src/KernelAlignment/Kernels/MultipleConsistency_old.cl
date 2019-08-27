#pragma once

#if defined SELECTIVITY_SUM
	#define SELECTIVITY_FUN(x,y) ((x)+(y))
#elif defined SELECTIVITY_MAX
	#define SELECTIVITY_FUN(x,y) max((x),(y))
#elif defined SELECTIVITY_MIN
	#define SELECTIVITY_FUN(x,y) min((x),(y))
#elif defined SELECTIVITY_AVG
	#define SELECTIVITY_FUN(x,y) (((x)+(y))/2)
#endif

#if defined FILTER_VERTICAL
	#define filter(a,b,x)			(((x) <= (a)) ? 1.0f : 0)
#elif defined FILTER_LINEAR
	#define filter(a,b,x)			((a) * (x) + (b))
#elif defined FILTER_HOMOGRAPHIC
	#define filter(a,b,x)			((1 - (x)) / ((a) * (x) + 1))
#elif defined FILTER_TRIANGLE_MIDPASS
	#define filter(a,b,x)			min((a)*(x), -(a)*(x)+(a))
#endif

/// <summary>
/// </summary>
/// <param name=lengths></param>
/// <param name=seqsWeights></param>
/// <param name=tasks></param>
/// <param name=sparseOffsets></param>
/// <param name=taskCount></param>
/// <param name=seqCount></param>
/// <param name=sparseMatrices></param>
/// <param name=posteriorMatrices></param>
/// <returns></returns>
__kernel __attribute__((reqd_work_group_size(STRIPE_COUNT * STRIPE_LENGTH, 1, 1)))
void MultipleConsistency_RelaxOld(
global	const	unsigned int* lengths,
global	const	unsigned int* sortingMap,
global	const	float* seqsWeights,
global	const	float* distances,
global	const	unsigned int* sparseOffsets,
global			struct RelaxationTask* tasks,
global			float* Sxy,						// output sparse matrix
global	const	float* Sxz,
global	const	float* Szy,
				struct RelaxationSchedule schedule,
local			struct cell_type* Sxz_elements)
{
	// make some aliases
	#define localId get_local_id(0)
	#define taskId get_group_id(0)
	#define stripeId (get_local_id(0) >> STRIPE_LENGTH_LOG2)
	
	int stripeOffset = get_local_id(0) % STRIPE_LENGTH; 
					// store stripeOffset and sparseWidth in the same variable	
	
	local index_type auxSizes[STRIPE_COUNT];
	local int x, y;
	local float wx_wy;
	local int length_x_inc;
	local int aux1;
	local int aux2;
	local float sumW;
	
	// distribution parameters
	local float a;
	local float b;
	
	if (localId == 0) {
		x = tasks[taskId].i; 
		y = tasks[taskId].j;
		length_x_inc = lengths[x] + 1;
		aux1 = sparseOffsets[x * schedule.numSeqs + y];
		sumW = *(Sxy + aux1); 
		
		a = distances[0 * schedule.numSeqs + 0];
		b = distances[1 * schedule.numSeqs + 1];

		// calculate adjusted selfweight sw = 1 + (sw - 1) * acceptedFraction
		wx_wy = 1.0f + (distances[schedule.numSeqs * schedule.numSeqs - 1] - 1.0f) * (float)tasks[taskId].acceptedCount / a;
		wx_wy *= seqsWeights[x] + seqsWeights[y];
	} 
	
	barrier(CLK_LOCAL_MEM_FENCE);

	Sxy += aux1;	// output sparse matrix

	// initialize local pointers and make some aliases
	local struct cell_type* Sxy_row = Sxz_elements + get_local_size(0);
	
	// add taks-specific offsets
	Sxz_elements += stripeId * STRIPE_LENGTH;
	Sxy_row += stripeId * schedule.sparseWidth;

	int end_ix = _ceilround(length_x_inc - 1, STRIPE_COUNT) + 1;
	int seed = tasks[taskId].seed;
	
	// contribution from all sequences other than x and y (uz - unsorted z index)
	for (int uz = schedule.sequenceBegin; uz < schedule.sequenceEnd; ++uz) {
		
		int z = sortingMap[uz];
		float w = SELECTIVITY_FUN(distances[x * schedule.numSeqs + z], distances[y * schedule.numSeqs + z]);
		
		seed = parkmiller(seed); // get next random number
		w = filter(a, b, w);
		w = ((float)seed) * RND_MAX_INV - w;  
		
		if (z == x || z == y || w >= 0) {
			continue;
		}
		
		w = seqsWeights[z] / (wx_wy);
		
		if (localId == 0) {
			aux1 = sparseOffsets[x * schedule.numSeqs + z];
			aux2 = sparseOffsets[z * schedule.numSeqs + y];
			sumW += w;
		}
		
		int length_z_inc = lengths[z] + 1;

		barrier(CLK_LOCAL_MEM_FENCE);
		Sxz += aux1;
		Szy += aux2;
		
		// iterate over Sxy rows
		for (int ix = 1 + stripeId; ix < end_ix; ix += STRIPE_COUNT) {
			
			// store Sxz sizes in auxiliary local array
			if (stripeOffset == 0) { 
				auxSizes[stripeId] = ix < length_x_inc ? sparse_rowSizes(Sxz)[ix] : 0; 
			}

			barrier(CLK_LOCAL_MEM_FENCE);
			int Sxy_size = ix < length_x_inc ? sparse_rowSizes(Sxy)[ix] : 0;
			int rowIndex = ix < length_x_inc ? sparse_rowIndices(Sxz, length_x_inc)[ix] : 0;
			int Sxz_size = auxSizes[stripeId]; 
			
			// get maximum Sxz length
			if (localId == 0) {
				index_type currentMax = 0;
				#pragma unroll STRIPE_COUNT
				for (int j = 0; j < STRIPE_COUNT; ++j) {
					currentMax = max(auxSizes[j], currentMax);
				} 
				auxSizes[0] = currentMax;
			}
			barrier(CLK_LOCAL_MEM_FENCE);
			
			// store locally Sxy rows being processed
			for (int y_lane = stripeOffset; y_lane < Sxy_size; y_lane += STRIPE_LENGTH) {
				Sxy_row[y_lane] = *sparse_cell(Sxy, length_x_inc, ix, y_lane);
			}
			
			--Sxy_size;
		
			barrier(CLK_LOCAL_MEM_FENCE);

			// process consecutive stripes of Sxz
			for (int z_lane = stripeOffset, end = _ceilround(auxSizes[0], STRIPE_LENGTH); z_lane < end; z_lane += STRIPE_LENGTH) {
				
				// store Sxz stripes locally
				if (z_lane < Sxz_size) { 
					Sxz_elements[stripeOffset] = sparse_data(Sxz, length_x_inc)[rowIndex + z_lane];
				}
				
				barrier(CLK_LOCAL_MEM_FENCE);
				
				// perform multiplication
				if (Sxy_size >= 0) 
				{
					int end_k = min(STRIPE_LENGTH, Sxz_size - z_lane + stripeOffset); // start from last index
					for (int k = 0; k < end_k ; ++k) {
									
						int Szy_size = sparse_rowSizes(Szy)[Sxz_elements[k].column];
						for (int y_lane = stripeOffset; y_lane < Szy_size; y_lane += STRIPE_LENGTH) {
			
							private struct cell_type Szy_elem = *sparse_cell(Szy, length_z_inc, Sxz_elements[k].column, y_lane);
							int lo = 0;
							int hi = Sxy_size;
								
							while (lo < hi) {
								int mid = hadd(lo, hi); // same as (lo + hi)>> 1;
								if (Sxy_row[mid].column >= Szy_elem.column) {
									hi = mid;
								} else {
									lo = mid + 1;
								}
							}

							if (Sxy_row[hi].column == Szy_elem.column) {	
								Sxy_row[hi].value += w * Sxz_elements[k].value * Szy_elem.value;
							}
						}
					}
				}
				
				barrier(CLK_LOCAL_MEM_FENCE);
			} // end stripe loop
			
			++Sxy_size;
			// copy Sxy rows back to global
			for (int y_lane = stripeOffset; y_lane < Sxy_size; y_lane += STRIPE_LENGTH) {
				sparse_cell(Sxy, length_x_inc, ix, y_lane)->value = Sxy_row[y_lane].value;
			}

			barrier(CLK_GLOBAL_MEM_FENCE);

		} // end row iteration
			
		Sxz -= aux1;
		Szy -= aux2;
	} // end relaxing sequence iteration

	if (localId == 0) {
		*Sxy = sumW;
		tasks[taskId].seed = seed;
	}

	#undef localId
	#undef taskId
	#undef stripeId
}

