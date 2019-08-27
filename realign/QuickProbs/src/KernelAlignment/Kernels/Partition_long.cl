#define GAP_OPEN		params->gapOpen
#define GAP_EXTEND		params->gapExt


#ifdef USE_FLOAT_PARTITION
	#define jagSize jagSize_float
	#define TERM_GAP_OPEN	1.0f
	#define TERM_GAP_EXTEND 1.0f
	#define PROB_LO 0.001f
	#define PROB_HI 1.0f
	
	typedef float partition_t;
#else
	#define jagSize jagSize_double
	#define TERM_GAP_OPEN	1.0
	#define TERM_GAP_EXTEND 1.0
	#define PROB_LO 0.001
	#define PROB_HI 1.0
	
	typedef double partition_t;
#endif

struct PartitionFunctionParams
{
	partition_t termGapOpen; 
	partition_t termGapExtend; 
	partition_t gapOpen; 
	partition_t gapExt;
	partition_t subMatrix[26 * 26];
};

/// <summary>
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=params></param>
/// <param name=aux></param>
/// <param name=forward></param>
/// <returns></returns>
__kernel 
void Partition_ComputeForward_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						struct PosteriorSchedule scheduling,
	PART_MEM	const	struct PartitionFunctionParams* params MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	local				partition_t* aux,
	global				partition_t* forward,
	global				partition_t* columns)
{
	
	#define j_offset get_local_id(0) 
	
	#define width (seq2Length + 1)
	#define height (seq1Length + 1)

	#define Zm forward
	#define Ze aux
	local partition_t* Zf = aux + get_local_size(0) * 1;
	local partition_t* Ze_Zf_Zm_prev = aux + get_local_size(0) * 2;

	#define Ze_column columns
	global partition_t* Ze_Zf_Zm_column = Ze_column + height;
	
	int j = j_offset;
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	partition_t Zm_up;

	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
		
		// column offset in global memory
		int ij_cur = jaggedIndex(j, -j_offset, height, jagSize);		// seek to -j'th row
		int ij_left = jaggedIndex(j - 1, -j_offset, height, jagSize);	
		
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		
		// initialization
		if (j <= seq2Length) { 
			Zf[j_offset] = 0; 
			Ze[j_offset] = (j > 0);
		} 
		
		if (j == 0) {
			Ze_Zf_Zm_prev[0] = 1; // next elements are not needed		
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		Zm_up = (j == 0);

		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset;	// calculate i coordinate starting from row 0
	
			partition_t Ze_Zf_Zm_diag;
			partition_t Ze_left;
			
			// utilise elements from Ze and Ze_Zf_Zm_prev arrays
			if (j > 0 && j < width && i >= 0 && i <= seq1Length) { 
				if (j_offset == 0) {
					Ze_Zf_Zm_diag = Ze_Zf_Zm_column[i];
					Ze_left = Ze_column[i];
				} else {
					Ze_Zf_Zm_diag = Ze_Zf_Zm_prev[j_offset - 1];
					Ze_left = Ze[j_offset - 1] ;
				}
			}
			barrier(CLK_LOCAL_MEM_FENCE);

			// update Ze_Zf_Zm_prev array
			if (j < width) {
				Ze_Zf_Zm_prev[j_offset] = Ze[j_offset] + Zf[j_offset] + Zm_up; 
			}

			// initializations for first column
			if (j == 0) { 
				if (i == 1) { 
					Zm_up = 0;
					Zf[0] = 1;  
				}
			} 

			barrier(CLK_LOCAL_MEM_FENCE);
		
			if (i >= 0 && i <= seq1Length && j <= seq2Length) // if indices are legal
			{
				// initialization of first row
				if (i > 0 && j > 0) {				
					int id = fast_mad(seq1[i], 26, seq2[j]);
					
					// update Zf 
					Zf[j_offset] = Zm_up * ((j == seq2Length) ? TERM_GAP_OPEN : GAP_OPEN) + 
						Zf[j_offset] * ((j == seq2Length) ? TERM_GAP_EXTEND : GAP_EXTEND);

					Zm_up = Ze_Zf_Zm_diag * params->subMatrix[id]; // now Zm_cur, future Zm_up
					Ze_Zf_Zm_diag = Zm[ij_left]; // use variable to temporarily store Zm_left
						
					// update Ze 
					Ze[j_offset] = Ze_Zf_Zm_diag * ((i == seq1Length) ? TERM_GAP_OPEN : GAP_OPEN) + 
						Ze_left * ((i == seq1Length) ? TERM_GAP_EXTEND: GAP_EXTEND);
				}
			
				Zm[ij_cur] = Zm_up;

				// last thread in a group stores column
				if (j_offset == get_local_size(0) - 1) {
					Ze_column[i] = Ze[j_offset];
					Ze_Zf_Zm_column[i] = Ze_Zf_Zm_prev[j_offset];
				}
			}

			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
			
			ij_cur += jagSize; // move to next row
			ij_left += jagSize;
		}	
	}

	//store the sum of zm zf ze (m,n)s in zm's 0,0 th position
	if (j - get_local_size(0) == seq2Length) {
		Zm[0] = Ze[j_offset] + Zf[j_offset] + Zm_up;
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);

	#undef height
	#undef width
	#undef Zm
	#undef j_offset
	#undef Ze_column
}


/// <summary>
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=params></param>
/// <param name=aux></param>
/// <param name=forward></param>
/// <param name=posterior></param>
/// <returns></returns>
__kernel 
void Partition_ComputeReverse_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						struct PosteriorSchedule scheduling,
	PART_MEM	const	struct PartitionFunctionParams* params MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	local				partition_t* aux,
	global	const		partition_t* forward,
	global				float* posterior,
	global				partition_t* columns)
{
	int j_offset = get_local_id(0);
	#define j_offset_rev (get_local_size(0) - get_local_id(0) - 1)
	#define width (seq2Length + 1)
	#define height (seq1Length + 1)
	
	++seq1; // + 1 because of '@' character at the beginning of seq.
	++seq2;

	#define Ze aux
	partition_t forward0 = forward[0];
	local partition_t* Zf = aux + 1 * get_local_size(0);
	local partition_t* Zm = aux + 2 * get_local_size(0);
	local partition_t* Ze_Zf_Zm_prev = aux + 3 * get_local_size(0);

	#define Ze_column columns
	global partition_t* Zm_column = Ze_column + height;
	global partition_t* Ze_Zf_Zm_column = Ze_column + 2 * height;

	int j = seq2Length - j_offset_rev;
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	
	barrier(CLK_GLOBAL_MEM_FENCE);

	for (int lane = 0; lane < lanesCount; ++lane, j -= get_local_size(0)) {

		int j_rev = seq2Length - j;
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		
		// calculate indices
		int ij_cur = jaggedIndex_float(j, seq1Length + j_offset_rev, height);					
		int ij_diag = jaggedIndex_float(j + 1, seq1Length + 1 + j_offset_rev, height);			
		int ij_forward_diag = jaggedIndex(j + 1, seq1Length + 1 + j_offset_rev, height, jagSize);

		if (j >= 0) {
			Zm[j_offset] = (j_rev == 0); // 0 ... 0 1
			Ze[j_offset] = (j_rev > 0); // 1 ... 1 0
			Zf[j_offset] = 0;
		}
	
		if (j_rev == 0)  {
			Ze_Zf_Zm_prev[j_offset] = 1;				
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		for (int iter = 0; iter < iterCount; ++iter) {
			int i = seq1Length + j_offset_rev - iter;
		
			partition_t Ze_Zf_Zm = 0;
			partition_t Zm_diag = 0;
			partition_t Ze_diag = 0;

			if (j_rev > 0 && i >= 0 && i <= seq1Length) { // except last column
				if (j_offset_rev == 0) {
					Ze_Zf_Zm = Ze_Zf_Zm_column[i];
					Ze_diag = Ze_column[i];
					Zm_diag = Zm_column[i];
				} else {
					Ze_Zf_Zm = Ze_Zf_Zm_prev[j_offset + 1];
					Zm_diag = Zm[j_offset + 1];
					Ze_diag = Ze[j_offset + 1];	
				}
			}

			barrier(CLK_LOCAL_MEM_FENCE);

			if (j_rev > 0 && j >= 0) {
				Ze_Zf_Zm_prev[j_offset] = Ze[j_offset] + Zf[j_offset] + Zm[j_offset];
			}
			else if (j_rev == 0)  {  // initialize last column
				if (i == seq1Length - 1) { 
					Zf[j_offset] = 1; 
					Zm[j_offset] = 0;
				}
			}
		
			barrier(CLK_LOCAL_MEM_FENCE);

			if (j >= 0 && i >= 0 && i <= seq1Length && j <= seq2Length) {
				
				#ifdef EXACT_REVERSE_PARTITION 
					if (i == 0 || j == 0) { posterior[ij_cur] = 0; } // initialization of first posterior row and column
				#endif

				if (i < seq1Length && j < seq2Length) // for all elements but last row and column
				{
					int id = fast_mad(seq1[i], 26, seq2[j]);

					// terminal gaps may be ignored here as last row and column are not used for posterior calculation
					Zf[j_offset] = Zm[j_offset] * GAP_OPEN + Zf[j_offset] * GAP_EXTEND; 
					Ze[j_offset] = Zm_diag * GAP_OPEN + Ze_diag * GAP_EXTEND;
				
					// use Zm_diag to store probability
					#define probability Zm_diag
					probability = (forward0 == 0) ? 0 : forward[ij_forward_diag] * Ze_Zf_Zm / forward0;
			
					Zm[j_offset] = Ze_Zf_Zm * params->subMatrix[id]; // multiply by score

					posterior[ij_diag] = (probability <= PROB_HI && probability >= PROB_LO) ? probability : 0;
					#undef probability
				}
				
				// first thread in a group stores column
				if (j_offset == 0) {
					Ze_column[i] = Ze[0];
					Zm_column[i] = Zm[0];
					Ze_Zf_Zm_column[i] = Ze_Zf_Zm_prev[0];
				}	
			}

			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

			ij_cur -= jagSize_float;
			ij_diag -= jagSize_float;
			ij_forward_diag -= jagSize;
		}	
	}		

	barrier(CLK_GLOBAL_MEM_FENCE);

	#undef height
	#undef width
	#undef j_offset_rev
	#undef Ze_column
}


/// <summary>
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=params></param>
/// <param name=aux></param>
/// <param name=posterior></param>
/// <returns></returns>
__kernel 
void Partition_ComputeAll_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						struct PosteriorSchedule scheduling,
	PART_MEM	const	struct PartitionFunctionParams* params MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	local				partition_t* auxiliaryRows, 
	global				partition_t* auxiliaryLayer,
	global				float* posterior,
	global				partition_t* columns)
{
	Partition_ComputeForward_long(seq1, seq2, seq1Length, seq2Length, scheduling, params, auxiliaryRows, auxiliaryLayer, columns);
	Partition_ComputeReverse_long(seq1, seq2, seq1Length, seq2Length, scheduling, params, auxiliaryRows, auxiliaryLayer, posterior, columns);
}
