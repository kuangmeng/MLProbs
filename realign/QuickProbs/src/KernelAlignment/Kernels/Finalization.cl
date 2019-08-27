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
void Finalization_CombineMatrices(
	global	const	float* matrix1, 
	global	const	float* matrix2,
					int seq1Length, 
					int seq2Length,
					float pid,
	local			float* localBuffer,
	global			float* out,
	global			float* extra) 
{
	#define j get_local_id(0)
	int columnSize = 0; // number of significant cells  

	local float* scoringRow = localBuffer;
	local int* columnSizes = (local int*)localBuffer;

//	local int* rowSizes = columnSizes + seq2Length + 1;
	// seq1 can be longer than the workgroup size
//	for (int q = j; q <= seq1Length; q += get_local_size(0)) {
//		rowSizes[q] = 0; 
//	}
	
	// column offset in global memory
	int ij = jaggedIndex_float(j, -j, seq1Length + 1);			// seek to -j'th row
	
	float diagScore = 0;
	float score = 0;
	
	for (int iter = 0, iterationsCount = seq1Length + seq2Length + 1; iter < iterationsCount; ++iter)
	{
		int i = iter - j;
		float leftScore = (j == 0 || j > seq2Length) ? 0 : scoringRow[j-1]; // calculate scoring matrix cells

		barrier(CLK_LOCAL_MEM_FENCE);

		if (i >= 0 && i <= seq1Length && j <= seq2Length) // if i and j are legal
		{	
			float a = 0;
			
			if (i > 0 && j > 0) {
				a = matrix1[ij];
				float b = matrix2[ij];
				a = sqrt(HMM_WEIGHT * a * a + (1 - HMM_WEIGHT) * b * b);					// store posterior cell in a
				
				score = _max3(a + diagScore, leftScore, score);		// store current score
				columnSize += (a >= POSTERIOR_CUTOFF) ? 1 : 0;
			}
			 
			out[ij] = a;
			scoringRow[j] = score;
			diagScore = leftScore;
		}
		
		// move to the next element
		ij += jagSize_float;
		
		barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
	}

	// store number of non zero cells in columns and rows
	if (j <= seq2Length)	{ columnSizes[j] = columnSize; }
	
	barrier(CLK_LOCAL_MEM_FENCE);

	// last thread sums up column sizes and saves stuff in extra output 
	if (j == seq2Length) {		
		for (int k = 1; k < seq2Length; ++k) { columnSize += columnSizes[k]; }
		
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
void Finalization_ComputeAll(
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
	global				float* posterior)
{
	#define j get_local_id(0)
	
	global float* probPosterior = auxiliaryThreeLayers;
	global float* funcPosterior = posterior;

	Partition_ComputeAll(seq1, seq2, seq1Length, seq2Length, scheduling, funcParams,
		localBuffer, auxiliaryThreeLayers, funcPosterior);
	
	Probabilistic_ComputeAll(seq1, seq2, seq1Length, seq2Length, scheduling, probParams,
		(local float*)localBuffer, probPosterior, probPosterior);

	#ifdef RETURN_GLOBAL
	if (j == 0) {
		probPosterior[0] = 0;
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
	#endif

	// calculate combined posterior matrix
	Finalization_CombineMatrices(
		probPosterior, 
		funcPosterior, 
		seq1Length, 
		seq2Length,
		pid,
		(local float*)localBuffer, 
		posterior, 
		probPosterior); // extra output will be placed at probPosterior
		
	// calculate distance on the basis of score
	if (j == 0) {
		posterior[0] = 1.0f - ((local float*)localBuffer)[0] / min(seq1Length, seq2Length);	
	}
	
	barrier(CLK_GLOBAL_MEM_FENCE);
	#undef j
}
