/// <summary>
/// Probabilistic parameters
/// </summary>
struct ProbabilisticParams
{
	float initial[NumMatrixTypes];						// holds the initial probabilities for each state
	float trans[NumMatrixTypes * NumMatrixTypes];        // holds all state-to-state transition probabilities
	float match[SymbolCountSmall * SymbolCountSmall];
	float insert[SymbolCountSmall * NumMatrixTypes];
};

#define PROBS_DIAG(j) probs->trans[fast_mad(j, depth, j)]

/// <summary>
/// Calculates forward matrix.
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=probabilistic></param>
/// <param name=levels></param>
/// <param name=forward></param>
/// <returns></returns>
__kernel 
void Probabilistic_ComputeForwardMatrix_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						struct PosteriorSchedule scheduling,
	HMM_MEM	const	struct ProbabilisticParams* probs MAX_SIZE(sizeof(struct ProbabilisticParams)),
	local				float* levels,
	global				float* forward,
	global				float* columns)		
{
	int height = seq1Length + 1;
	int width = seq2Length + 1;
	int j_offset = get_local_id(0);
	int j = j_offset;
	#define layerOffset ((int)get_local_size(0))

	for (int q = j_offset; q < depth * height; q += get_local_size(0)) {
		columns[q] = LOG_ZERO;
	}

	barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
	
	private float levels_up[NumInsertStates];
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	float fwd_up;
	float levels_diag;

	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
		
		// calculate indices
		int c2 = j <= seq2Length ? seq2[j] : -1;
		int ij_cur = jaggedIndex_float(j, -j_offset, height);		// 0 - current					
		int ij_left = jaggedIndex_float(j - 1, -j_offset, height);	// 2 - left
	
		fwd_up = LOG_ZERO;
		levels_diag = LOG_ZERO;
		for (int k = 0; k < NumInsertStates; ++k) {
			levels_up[k] = LOG_ZERO;
		}
		
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		
		for (int iter = 0; iter < iterCount; ++iter)
		{
			private float levels_left[NumInsertStates];
			float fwd_cur;
			float fwd_left;
		
			int c1;
			int i = iter - j_offset;						// calculate i coordinate	
			bool legal = i >= 0 && i < height && j < width;  // if i and j are legal
		
			if (legal) {
				c1 = seq1[i];
				int ixj = i * j;

			fwd_cur = (ixj == 0) ? 
					LOG_ZERO :
					(ixj == 1 ? probs->initial[0] : levels_diag) + probs->match[fast_mad(c1, SymbolCountSmall, c2)];
			fwd_left = (j > 0) ? forward[ij_left] : LOG_ZERO;
			forward[ij_cur] = fwd_cur;

			// store elements in left cells, init elements in up_cells from 0th column
			for (int k = 0; k < NumInsertStates; ++k) {
				levels_up[k] = (i == 1 && j == 0) ?
					probs->initial[fast_mad(2,k,1)] + probs->insert[fast_mad(c1, depth, k)] :
					levels_up[k]; 
				
				levels_left[k] = (i == 0 && j == 1) ? 
					probs->initial[fast_mad(2,k,2)] + probs->insert[fast_mad(c2, depth, k)] :
					((j > 0) ? 
						(j_offset == 0) ? 
							columns[k * height + i] : 
							levels[fast_mad(layerOffset, k, j_offset - 1)] : LOG_ZERO);
				}

				levels_diag = (j_offset == 0) ? columns[2 * height + i] : levels[fast_mad(layerOffset, 2, j_offset - 1)];
			}

			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE); // synchronize after utilising levels elements
		
			bool legalEx = legal && (i > 1 || j > 1); // alter 'legal 'flag
			if (legalEx) {
				int q = 1;
				for (int k = 0; k < NumInsertStates; ++k) {
					levels_up[k] = (i > 0) ?
						probs->insert[fast_mad(c1, depth, k)] + logOfSum(fwd_up + probs->trans[q], levels_up[k] + PROBS_DIAG(q)) :
						levels_up[k];
					++q;
					levels_left[k] = (j > 0) ? 
						probs->insert[fast_mad(c2, depth, k)] + logOfSum(fwd_left + probs->trans[q], levels_left[k] + PROBS_DIAG(q)) :
						levels_left[k];
					++q;
				}
			}
		
			// update all local memory rows
			if (legal) {
				fwd_left = fwd_cur + probs->trans[0]; // use fwd_left as auxiliary 
				for (int k = 0; k < NumInsertStates; ++k) {
					logOfSum_private(&fwd_left, levels_up[k] + probs->trans[fast_mul(fast_mad(2,k,1), depth)]); // current cell in up layers
					logOfSum_private(&fwd_left, levels_left[k] + probs->trans[fast_mul(fast_mad(2,k,2), depth)]); // current cell in left layers
					levels[fast_mad(layerOffset, k, j_offset)] = levels_left[k];
				}
				levels[fast_mad(layerOffset, 2, j_offset)] = fwd_left;
			}

			fwd_up = fwd_cur; 
			
			// compute total probability from right-bottom element of forward matrix
			// (corresponding fragment of backward matrix is initialized with known values
			// so we can use them directly instead of backward matrix itself)
			if (i == seq1Length && j == seq2Length) {
				fwd_left = fwd_cur + probs->initial[0]; // use fwd_left for storing total probability		
				for (int k = 0; k < NumInsertStates; ++k) {
					logOfSum_private(&fwd_left, levels_up[k] + probs->initial[fast_mad(2,k,1)]); // up
					logOfSum_private(&fwd_left, levels_left[k] + probs->initial[fast_mad(2,k,2)]); // components from other layers
				}
			
				levels[0] = fwd_left;
			
				#ifdef RETURN_GLOBAL
					forward[0] = fwd_left;
				#endif
			}

			barrier(CLK_LOCAL_MEM_FENCE);

			if (i >= 0 && i < height && (j_offset == get_local_size(0) - 1)) {
				for (int k = 0; k < NumInsertStates + 1; ++k) { // as there are three layers
					columns[k * height + i] = levels[fast_mad(layerOffset, k, j_offset)];
				}
			}
		
			ij_cur += jagSize_float; // move to next row
			ij_left += jagSize_float; // move to next row
		
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE); // synchronize after every element

		}	// end vertical loop
	}

	barrier(CLK_LOCAL_MEM_FENCE);
	#undef layerOffset
}


/// <summary>
/// Calculates backward matrix.
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=probabilistic></param>
/// <param name=levels></param>
/// <param name=backward></param>
/// <returns></returns>
__kernel 
void Probabilistic_ComputeBackwardMatrix_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						struct PosteriorSchedule scheduling,
	HMM_MEM	const	struct ProbabilisticParams* probs MAX_SIZE(sizeof(struct ProbabilisticParams)),
	local				float* levels,
	global				float* backward,
	global				float* columns)		
{
	int height = seq1Length + 1;
	int width = seq2Length + 1;
	int j_offset = get_local_id(0);
	#define j_offset_rev (get_local_size(0) - get_local_id(0) - 1)
	int j = seq2Length - j_offset_rev;
	#define layerOffset ((int)get_local_size(0))

	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	float total = 0;

	for (int lane = 0; lane < lanesCount; ++lane, j -= get_local_size(0)) {
		int j_rev = seq2Length - j;
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
		
		// remember offset for each index combination
		int ij_cur = jaggedIndex_float(j, seq1Length + j_offset_rev, height);						
		int ij_diag = jaggedIndex_float(j + 1, seq1Length + j_offset_rev, height) + jagSize_float;
	
		// initialization condition

		private float levels_down[NumInsertStates];
		private float levels_right[NumInsertStates];
		for (int k = 0; k < NumInsertStates; ++k) {
			levels_down[k] = j_rev == 0 ? probs->initial[fast_mad(2,k,1)] : LOG_ZERO;
			levels_right[k] = j_rev == 0 ? probs->initial[fast_mad(2,k,2)] : LOG_ZERO; // this will be used for current and right
		}

		for (int iter = 0; iter < iterCount; ++iter) {
			
			int i = seq1Length + j_offset_rev - iter;
			bool legal = i >= 0 && i <= seq1Length && j >= 0 && j < width;
			if (legal) {
				int c1 = seq1[i + 1]; // + 1 as in original code (assume seq1 buffer to be longer)
				int c2 = seq2[j + 1];
		
				float bwd_cur = (i == seq1Length && j == seq2Length) ? probs->initial[0] : LOG_ZERO;
				float bwd_diag_ex = LOG_ZERO;
			
				if (i < seq1Length && j < seq2Length) {
					bwd_diag_ex = backward[ij_diag] + probs->match[fast_mad(c1, SymbolCountSmall, c2)];
					bwd_cur = bwd_diag_ex + probs->trans[0];
				}

				// calculate total probability from left-top fragment of backward matrix
				// (corresponding fragment of forward matrix is initialized with known values 
				// so we can use them directly instead of forward matrix itself)
				if (i == 0 && j == 0) {
					total = probs->initial[0] + bwd_diag_ex;
					for (int k = 0; k < NumInsertStates; ++k) {
						logOfSum_private(&total, probs->initial[fast_mad(2,k,1)] + probs->insert[fast_mad(c1, depth, k)] + levels_down[k]);
						logOfSum_private(&total, probs->initial[fast_mad(2,k,2)] + probs->insert[fast_mad(c2, depth, k)] + 
							((j_offset_rev == 0) ? columns[k * height + i] : levels[fast_mad(layerOffset, k, j_offset + 1)]));							
					}
				}

				if (i < seq1Length) {
					for (int k = 0; k < NumInsertStates; ++k) {
						int q = fast_mad(2,k,1);
						levels_right[k] = LOG_ZERO;
			
						float aux = levels_down[k] + probs->insert[fast_mad(c1, depth, k)]; 
						levels_down[k] = (j < seq2Length) ? bwd_diag_ex + probs->trans[fast_mul(q, depth)] : LOG_ZERO;
						logOfSum_private(&bwd_cur, aux + probs->trans[q]);
						logOfSum_private(&levels_down[k], aux + probs->trans[fast_mad(q, depth, q)]);
					}
				}
			
				if (j < seq2Length) {
					for (int k = 0; k < NumInsertStates; ++k) {
						int q = fast_mad(2,k,2);
						float aux = 
							((j_offset_rev == 0) ? columns[k * height + i] : levels[fast_mad(layerOffset, k, j_offset + 1)])
							 + probs->insert[fast_mad(c2, depth, k)];
						
						levels_right[k] = (i < seq1Length) ? bwd_diag_ex + probs->trans[fast_mul(q, depth)] : LOG_ZERO;
						logOfSum_private(&bwd_cur, aux + probs->trans[q]);
						logOfSum_private(&levels_right[k], aux + probs->trans[fast_mad(q, depth, q)]);
					}
				}

				backward[ij_cur] = bwd_cur;
			}
		
			barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);

			// first thread in a group stores column
			if (legal) {
				for (int k = 0; k < NumInsertStates; ++k) {
					if (j_offset == 0 ) {
						columns[k * height + i] = levels_right[k];
					} else {
						levels[fast_mad(layerOffset, k, j_offset)] = levels_right[k]; 		// store private cells
					}
				}	
			}
				
			ij_cur -= jagSize_float;
			ij_diag -= jagSize_float;
		
			barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);
		} // end vertical loop
	}

	if (total != 0) {
		levels[0] = total;
		#ifdef RETURN_GLOBAL
			backward[jaggedIndex(seq2Length, seq1Length, height, jagSize_float)] = total;
		#endif
	}

	barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

	#undef levelWidth
}


/// <summary>
/// Calculates total probability on the basis of forward and backward matrices.
/// </summary>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=totalProb></param>
/// <param name=scheduling></param>
/// <param name=forward></param>
/// <param name=backward></param>
/// <param name=posterior></param>
/// <returns></returns>
__kernel 
void Probabilistic_ComputePosteriorMatrix_long(       
					int seq1Length, 
					int seq2Length,
					float totalProb,
					struct PosteriorSchedule scheduling,
	global	const	float* forward, 
	global	const	float* backward,
	global			float* posterior)
{
	int j_offset = get_local_id(0);
	int lanesCount = _ceildiv(seq2Length + 1, get_local_size(0));
	totalProb = (totalProb == 0) ? 1.00f : totalProb;

	int j = j_offset;
	for (int lane = 0; lane < lanesCount; ++lane, j += get_local_size(0)) {
		
		int ij = jaggedIndex_float(j, -j_offset, seq1Length + 1);
		int iterCount = seq1Length + min(get_local_size(0), (uint)(seq2Length + 1) - lane * get_local_size(0));
	
		for (int iter = 0; iter < iterCount; ++iter) {
			int i = iter - j_offset;
			if (i >= 0 && i <= seq1Length && j <= seq2Length) {
				posterior[ij] = fast_exp(min(LOG_ONE, forward[ij] + backward[ij] - totalProb));
			}
			ij += jagSize_float;
		}
	
		barrier(CLK_GLOBAL_MEM_FENCE);
	}
}

					
/// <summary>
/// Calculates probabilistic posterior matrix using sequences directly.
/// </summary>
/// <param name=seq1></param>
/// <param name=seq2></param>
/// <param name=seq1Length></param>
/// <param name=seq2Length></param>
/// <param name=scheduling></param>
/// <param name=probabilistic></param>
/// <param name=buffer_2layers></param>
/// <param name=posterior></param>
/// <returns></returns>
__kernel 
void Probabilistic_ComputeAll_long(
	global		const	char* seq1,
	global		const	char* seq2,
						int seq1Length,
						int seq2Length,
						struct PosteriorSchedule scheduling,
	HMM_MEM	const	struct ProbabilisticParams* probabilistic MAX_SIZE(sizeof(struct ProbabilisticParams)),
	local				float* auxiliaryRows,
	global				float* auxiliaryTwoLayers,
	global				float* posterior,
	global				float* columns)
{
	#define j get_local_id(0)
	int layerSize = jaggedSize(seq2Length + 1, seq1Length + 1, jagSize_float);
	
	// place forward matrix at layer 0, backward matrix at layer 1, treat posterior as auxiliary buffer
	#define forward	(auxiliaryTwoLayers)
	#define backward (auxiliaryTwoLayers + layerSize)
	
	Probabilistic_ComputeForwardMatrix_long(
		seq1, seq2, seq1Length, seq2Length, scheduling, probabilistic, auxiliaryRows, forward, columns);
	float totalProb = auxiliaryRows[0]; // get total probability from forward algorithm
	barrier(CLK_LOCAL_MEM_FENCE);

	Probabilistic_ComputeBackwardMatrix_long(
		seq1, seq2, seq1Length, seq2Length, scheduling, probabilistic, auxiliaryRows, backward, columns);
	totalProb = (totalProb + auxiliaryRows[0]) / 2; // get total probability from backward algorithm and combine it with forward-one
	barrier(CLK_LOCAL_MEM_FENCE);

	#ifdef RETURN_GLOBAL
	// restore original values
	if (j == 0) {
		int lastIndex = jaggedIndex_float(seq2Length, seq1Length, seq1Length + 1);
		forward[0] = LOG_ZERO;
		backward[lastIndex] = probabilistic->initial[0];
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
	#endif

	// compute posterior probability
	Probabilistic_ComputePosteriorMatrix_long(seq1Length, seq2Length, totalProb, scheduling, forward, backward, posterior);

	// This is just to return total probability. Remember to set it to 0 in host! 
	#ifdef RETURN_GLOBAL
	if (j == 0) {
		posterior[0] = totalProb;
	}
	barrier(CLK_GLOBAL_MEM_FENCE);
	#endif

	#undef j
	#undef forward
	#undef backward
}
