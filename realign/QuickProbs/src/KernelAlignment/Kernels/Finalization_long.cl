/// <summary>
/// Calculates total probability on the basis of forward and backward matrices.
/// </summary>
/// <param name=matrix1></param>
/// <param name=matrix2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=buffer_row></param>
/// <param name=out></param>
/// <param name=score></param>
/// <returns></returns>
__kernel 
void Finalization_CombineMatrices_long(
	global	const	float* matrix1, 
	global	const	float* matrix2,
					int seq1Length, 
					int seq2Length,
					float pid,
	local			float* localBuffer,
	global			float* out,
	global			float* extra,
	global			float* scoringColumn) 
{
	#define j_offset get_local_id(0) 
	int columnSize = 0; // number of significant cells  

	local float* scoringRow = localBuffer;
	local int* columnSizes = (local int*)localBuffer;

	// reset scoring column
	for (int i = get_local_id(0); i <= seq1Length; i += get_local_size(0)) {
		scoringColumn[i] = 0;
	}
	 
	barrier(CLK_GLOBAL_MEM_FENCE);
	
	float score = 0;
	
	int j = j_offset;
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
		
		float diagScore = 0;
		score = 0;

		// column offset in global memory
		int ij = jaggedIndex_float(j, -j_offset, seq1Length + 1);			// seek to -j'th row
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		

		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset;

			if (i >= 0 && i <= seq1Length && j <= seq2Length) { 
				// if i and j are legal
				float leftScore = (j_offset == 0) ? scoringColumn[i] : scoringRow[j_offset - 1];
				
				float a = 0;
			
				if (i > 0 && j > 0) {
					a = matrix1[ij];
					float b = matrix2[ij];
					a = sqrt(HMM_WEIGHT * a * a + (1 - HMM_WEIGHT) * b * b); // store posterior cell in a
				
					score = _max3(a + diagScore, leftScore, score);		// store current score
					columnSize += (a >= POSTERIOR_CUTOFF) ? 1 : 0;
				}
			 
				out[ij] = a;
				scoringRow[j_offset] = score;
				diagScore = leftScore;

				// last thread in a group stores column
				if (j_offset == get_local_size(0) - 1) {
					scoringColumn[i] = score;
				}
			}

			// move to the next element
			ij += jagSize_float;
		
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		}
	}

	// store number of non zero cells in columns and rows
	columnSizes[get_local_id(0)] = columnSize;
	barrier(CLK_LOCAL_MEM_FENCE);

	// last saves stuff in extra output 
	if (j - get_local_size(0) == seq2Length) {		
		
		columnSize = 0;
		for (int k = 0; k < get_local_size(0); ++k) { 
			columnSize += columnSizes[k]; 
		}
		
		localBuffer[0] = score;
		localBuffer[1] = columnSize;
		
		#ifdef RETURN_GLOBAL
		if (extra != 0) {
			extra[0] = score;
			extra[1] = columnSize; 
		}
		#endif
	} 
	
	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
}

/// <summary>
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=probParams></param>
/// <param name=funcParams></param>
/// <param name=buffer_3layers></param>
/// <param name=posterior></param>
/// <returns></returns>
__kernel 
void Finalization_ComputeAll_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						float pid,
						struct PosteriorSchedule scheduling,
	HMM_MEM		const	struct ProbabilisticParams* probParams,
	PART_MEM	const	struct PartitionFunctionParams* funcParams,
	local				partition_t* localBuffer,
	global				partition_t* auxiliaryThreeLayers,
	global				float* posterior,
	global				float* columns)
{
	#define j get_local_id(0)
	
	global float* probPosterior = auxiliaryThreeLayers;
	global float* funcPosterior = posterior;

	Partition_ComputeAll_long(seq1, seq2, seq1Length, seq2Length, scheduling, funcParams,
		localBuffer, auxiliaryThreeLayers, funcPosterior, columns);
	
	Probabilistic_ComputeAll_long(seq1, seq2, seq1Length, seq2Length, scheduling, probParams,
		(local float*)localBuffer, probPosterior, probPosterior, columns);

	#ifdef RETURN_GLOBAL
	if (j == 0) {
		probPosterior[0] = 0;
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
	#endif

	// calculate combined posterior matrix
	Finalization_CombineMatrices_long(
		probPosterior, 
		funcPosterior, 
		seq1Length, 
		seq2Length,
		pid,
		(local float*)localBuffer, 
		posterior, 
		probPosterior, // extra output will be placed at probPosterior
		columns); 
		
	// calculate distance on the basis of score
	if (j == 0) {
		posterior[0] = 1.0f - ((local float*)localBuffer)[0] / min(seq1Length, seq2Length);	
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
	#undef j
}

