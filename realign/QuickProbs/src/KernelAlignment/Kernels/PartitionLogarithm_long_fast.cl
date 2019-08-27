#define GAP_OPEN		params->gapOpen
#define GAP_EXTEND		params->gapExt

#define TERM_GAP_OPEN	LOG_ONE
#define TERM_GAP_EXTEND LOG_ONE
#define PROB_LO 0.001f
#define PROB_HI 1.0f
typedef float partition_t;

struct PartitionFunctionParams
{
	float termGapOpen; 
	float termGapExtend; 
	float gapOpen; 
	float gapExt;
	float subMatrix[26 * 26];
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
	local				float* aux,
	global				float* forward,
	global				float* columns)
{
	
	#define j_offset get_local_id(0) 
	#define width (seq2Length + 1)
	#define height (seq1Length + 1)

	#define Zm forward
	#define Ze aux
	local float* Ze_Zf_Zm_prev = aux + get_local_size(0);

	#define Ze_column columns
	global float* Ze_Zf_Zm_column = Ze_column + height;
	
	int j = j_offset;
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
		
		// column offset in global memory
		int ij_cur = jaggedIndex(j, -j_offset, height, JAG_FLOAT);		// seek to -j'th row
		int ij_left = jaggedIndex(j - 1, -j_offset, height, JAG_FLOAT);	
		
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		
		// initialization
		if (j <= seq2Length) { 
			Ze[j_offset] = (j > 0) ? LOG_ONE : LOG_ZERO;
		} 

		barrier(CLK_LOCAL_MEM_FENCE);

		float Zf_up = LOG_ZERO;
		float Zm_up = (j == 0) ? LOG_ONE : LOG_ZERO;

		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset;	// calculate i coordinate starting from row 0
	
			float Ze_Zf_Zm_diag;
			float Ze_left;
			bool legal = i >= 0 && i <= seq1Length && j <= seq2Length;

			// utilise elements from Ze and Ze_Zf_Zm_prev arrays
			if (legal && j > 0) { 
				if (j_offset == 0) {
					Ze_left = Ze_column[i];
				} else {
					Ze_left = Ze[j_offset - 1] ;
				}
			}
			
			// initializations for first column
			if (j == 0 && i == 1) { 
				Zm_up = LOG_ZERO;
				Zf_up = LOG_ONE;  
			} 

			barrier(CLK_LOCAL_MEM_FENCE);
		
			if (legal) { // if indices are legal
				// initialization of first row
				if (i > 0 && j > 0) {				
					int id = fast_mad(seq1[i], 26, seq2[j]);
					
					// update Zf 
					float t = Zm_up + ((j == seq2Length) ? TERM_GAP_OPEN : GAP_OPEN);
					logOfSum_private(&t, Zf_up + ((j == seq2Length) ? TERM_GAP_EXTEND : GAP_EXTEND));
					Zf_up = t; 
				
					Zm_up = Ze_Zf_Zm_diag + params->subMatrix[id]; // current Zm_s0, future Zm_up
					Ze_Zf_Zm_diag = Zm[ij_left]; // use variable to temporarily store Zm_left
						
					// update Ze 
					t = Ze_Zf_Zm_diag + ((i == seq1Length) ? TERM_GAP_OPEN : GAP_OPEN);
					logOfSum_private(&t, Ze_left + ((i == seq1Length) ? TERM_GAP_EXTEND: GAP_EXTEND));
					Ze[j_offset] = t;
				}
			
				Zm[ij_cur] = Zm_up;
				Ze_Zf_Zm_diag = (j > 0 && j_offset == 0) ? Ze_Zf_Zm_column[i] : Ze_Zf_Zm_prev[j_offset - 1]; // load prev elements
			}
				
			barrier(CLK_LOCAL_MEM_FENCE);

			if (legal) {
				// update Ze_Zf_Zm_prev array
				float t; 
				t = Ze[j_offset];
				logOfSum_private(&t, Zf_up);
				logOfSum_private(&t, Zm_up);
				Ze_Zf_Zm_prev[j_offset] = t;
				
				// last thread in a group stores column
				if (j_offset == get_local_size(0) - 1) {
					Ze_column[i] = Ze[j_offset];
					Ze_Zf_Zm_column[i] = t;
				} 
			}

			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
			
			ij_cur += JAG_FLOAT; // move to next row
			ij_left += JAG_FLOAT;
		}	
	}

	//store the sum of zm zf ze (m,n)s in zm's 0,0 th position
	if (j - get_local_size(0) == seq2Length) {
		Zm[0] = Ze_Zf_Zm_prev[j_offset];
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
	local				float* aux,
	global	const		float* forward,
	global				float* posterior,
	global				float* columns)
{
	int j_offset = get_local_id(0);
	#define j_offset_rev (get_local_size(0) - get_local_id(0) - 1)
	#define width (seq2Length + 1)
	#define height (seq1Length + 1)
	
	++seq1; // + 1 because of '@' character at the beginning of seq.
	++seq2;

	#define Ze aux
	float forward0 = forward[0];
	local float* Zm = aux + 1 * get_local_size(0);
	local float* Ze_Zf_Zm_prev = aux + 2 * get_local_size(0);

	#define Ze_column columns
	global float* Zm_column = Ze_column + height;
	global float* Ze_Zf_Zm_column = Ze_column + 2 * height;

	int j = seq2Length - j_offset_rev;
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	
	barrier(CLK_GLOBAL_MEM_FENCE);

	for (int lane = 0; lane < lanesCount; ++lane, j -= get_local_size(0)) {

		int j_rev = seq2Length - j;
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		
		// calculate indices
		int ij_cur = jaggedIndex(j, seq1Length + j_offset_rev, height, JAG_FLOAT);					
		int ij_diag = jaggedIndex(j + 1, seq1Length + 1 + j_offset_rev, height, JAG_FLOAT);			
		int ij_forward_diag = jaggedIndex(j + 1, seq1Length + 1 + j_offset_rev, height, JAG_FLOAT);

		float Zf = LOG_ZERO;

		if (j >= 0) {
			Zm[j_offset] = (j_rev == 0) ? LOG_ONE : LOG_ZERO; // 0 ... 0 1
			Ze[j_offset] = (j_rev > 0) ? LOG_ONE : LOG_ZERO; // 1 ... 1 0
		}
	
		if (j_rev == 0)  {
			Ze_Zf_Zm_prev[j_offset] = LOG_ONE;				
		}

		barrier(CLK_LOCAL_MEM_FENCE);

		float Ze_Zf_Zm = LOG_ZERO;

		for (int iter = 0; iter < iterCount; ++iter) {
			int i = seq1Length + j_offset_rev - iter;
			
			float Zm_diag = LOG_ZERO;
			float Ze_diag = LOG_ZERO;
			bool legal = j >= 0 && i >= 0 && i <= seq1Length && j <= seq2Length;

			if (legal && j_rev > 0) { // except last column
				if (j_offset_rev == 0) {
					Ze_diag = Ze_column[i];
					Zm_diag = Zm_column[i];
				} else {
					Zm_diag = Zm[j_offset + 1];
					Ze_diag = Ze[j_offset + 1];	
				}
			}

			if (j_rev == 0 && i == seq1Length - 1)  {  // initialize last column
				Zf = LOG_ONE; 
				Zm[j_offset] = LOG_ZERO;
			}
		
			barrier(CLK_LOCAL_MEM_FENCE);

			if (legal) {
				
				#ifdef EXACT_REVERSE_PARTITION 
					if (i == 0 || j == 0) { posterior[ij_cur] = 0; } // initialization of first posterior row and column
				#endif

				if (i < seq1Length && j < seq2Length) // for all elements but last row and column
				{
					int id = fast_mad(seq1[i], 26, seq2[j]);

					// terminal gaps may be ignored here as last row and column are not used for posterior calculation
					float t = Zm[j_offset] + GAP_OPEN;
					logOfSum_private(&t, Zf + GAP_EXTEND); 
					Zf = t;

					t = Zm_diag + GAP_OPEN;
					logOfSum_private(&t, Ze_diag + GAP_EXTEND);
					Ze[j_offset] = t; 
				
					// use t to store probability
					t = (forward0 == LOG_ZERO) ? 0 : fast_exp(forward[ij_forward_diag] + Ze_Zf_Zm - forward0);
					Zm[j_offset] = Ze_Zf_Zm + params->subMatrix[id]; // multiply by score
					posterior[ij_diag] = (t <= PROB_HI && t >= PROB_LO) ? t : 0;
				}

				Ze_Zf_Zm = (j_offset_rev == 0) ?  (j_rev == 0 ? LOG_ZERO : Ze_Zf_Zm_column[i]) : Ze_Zf_Zm_prev[j_offset + 1];
			}
				
			barrier(CLK_LOCAL_MEM_FENCE);

			if (legal) {
				float t;
				if (j_rev > 0 && j >= 0) {
					t = Ze[j_offset];
					logOfSum_private(&t, Zf);
					logOfSum_private(&t, Zm[j_offset]);
					Ze_Zf_Zm_prev[j_offset] = t;
				}

				// first thread in a group stores column
				if (j_offset == 0) {
					Ze_column[i] = Ze[0];
					Zm_column[i] = Zm[0];
					Ze_Zf_Zm_column[i] = t;
				}	
			}

			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

			ij_cur -= JAG_FLOAT;
			ij_diag -= JAG_FLOAT;
			ij_forward_diag -= JAG_FLOAT;
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
	local				float* auxiliaryRows, 
	global				float* auxiliaryLayer,
	global				float* posterior,
	global				float* columns)
{
	Partition_ComputeForward_long(seq1, seq2, seq1Length, seq2Length, scheduling, params, auxiliaryRows, auxiliaryLayer, columns);
	Partition_ComputeReverse_long(seq1, seq2, seq1Length, seq2Length, scheduling, params, auxiliaryRows, auxiliaryLayer, posterior, columns);
}
