#pragma once

#ifndef LONG_KERNELS

__kernel void MultipleProbabilistic_ComputeAll(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localMemory)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PROBABILISTIC_PARAMS
		// copy scoring matrices from global to local memory
		local struct ProbabilisticParams localParams;
		local struct ProbabilisticParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct ProbabilisticParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Probabilistic_ComputeAll(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij,
		auxiliaryMatrices + tasks[taskId].offset_Aij);
}


__kernel void MultipleProbabilistic_ComputeForward(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localMemory)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PROBABILISTIC_PARAMS
		// copy scoring matrices from global to local memory
		local struct ProbabilisticParams localParams;
		local struct ProbabilisticParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct ProbabilisticParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Probabilistic_ComputeForwardMatrix(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij);
}


__kernel void MultipleProbabilistic_ComputeBackward(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localMemory)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	int length_i = lengths[task_i];
	int length_j = lengths[task_j];

	int layerSize = jaggedSize(length_j + 1, length_i + 1, jagSize_float);

	#ifdef COPY_PROBABILISTIC_PARAMS
		// copy scoring matrices from global to local memory
		local struct ProbabilisticParams localParams;
		local struct ProbabilisticParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct ProbabilisticParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Probabilistic_ComputeBackwardMatrix(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		length_i,
		length_j,
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij + layerSize) ;
}

__kernel void MultipleProbabilistic_Combine(
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	int length_i = lengths[task_i];
	int length_j = lengths[task_j];

	int layerSize = jaggedSize(length_j + 1, length_i + 1, jagSize_float);
	
	global float* forward	= auxiliaryMatrices + tasks[taskId].offset_Aij;
	global float* backward = auxiliaryMatrices + tasks[taskId].offset_Aij + layerSize;

	local float localTotal;

	// restore original values
	if (get_local_id(0) == 0) {
		int lastIndex = jaggedIndex_float(length_j, length_i, length_i + 1);
		localTotal = (forward[0] + backward[lastIndex]) / 2;
		forward[0] = LOG_ZERO;
		backward[lastIndex] = globalParams->initial[0];
	}

	barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

	float total = localTotal;
	
	Probabilistic_ComputePosteriorMatrix(
		length_i,
		length_j,
		total,
		scheduling,
		forward,
		backward,
		auxiliaryMatrices + tasks[taskId].offset_Aij);
}


#else

__kernel void MultipleProbabilistic_ComputeAll_long(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localMemory,
	global				float*			columns,
				const	int				columnsOffset)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PROBABILISTIC_PARAMS
		// copy scoring matrices from global to local memory
		local struct ProbabilisticParams localParams;
		local struct ProbabilisticParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct ProbabilisticParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Probabilistic_ComputeAll_long(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij,
		auxiliaryMatrices + tasks[taskId].offset_Aij,
		columns + taskId * columnsOffset);
}


__kernel void MultipleProbabilistic_ComputeForward_long(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localMemory,
	global				float*			columns,
				const	int				columnsOffset)
				
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PROBABILISTIC_PARAMS
		// copy scoring matrices from global to local memory
		local struct ProbabilisticParams localParams;
		local struct ProbabilisticParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct ProbabilisticParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Probabilistic_ComputeForwardMatrix_long(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij,
		columns + taskId * columnsOffset);
}

__kernel void MultipleProbabilistic_ComputeBackward_long(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localMemory,
	global				float*			columns,
				const	int				columnsOffset)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	int length_i = lengths[task_i];
	int length_j = lengths[task_j];

	int layerSize = jaggedSize(length_j + 1, length_i + 1, jagSize_float);

	#ifdef COPY_PROBABILISTIC_PARAMS
		// copy scoring matrices from global to local memory
		local struct ProbabilisticParams localParams;
		local struct ProbabilisticParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct ProbabilisticParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Probabilistic_ComputeBackwardMatrix_long(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		length_i,
		length_j,
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij + layerSize,
		columns + taskId * columnsOffset) ;
}

__kernel void MultipleProbabilistic_Combine_long(
	global		const	unsigned int*	lengths,
	global		const	struct			ProbabilisticParams* globalParams MAX_SIZE(sizeof(struct ProbabilisticParams)),
	global				float*			auxiliaryMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	int length_i = lengths[task_i];
	int length_j = lengths[task_j];

	int layerSize = jaggedSize(length_j + 1, length_i + 1, jagSize_float);
	
	global float* forward	= auxiliaryMatrices + tasks[taskId].offset_Aij;
	global float* backward = auxiliaryMatrices + tasks[taskId].offset_Aij + layerSize;

	local float localTotal;

	// restore original values
	if (get_local_id(0) == 0) {
		int lastIndex = jaggedIndex_float(length_j, length_i, length_i + 1);
		localTotal = (forward[0] + backward[lastIndex]) / 2;
		forward[0] = LOG_ZERO;
		backward[lastIndex] = globalParams->initial[0];
	}

	barrier(CLK_GLOBAL_MEM_FENCE | CLK_LOCAL_MEM_FENCE);

	float total = localTotal;
	
	Probabilistic_ComputePosteriorMatrix_long(
		length_i,
		length_j,
		total,
		scheduling,
		forward,
		backward,
		forward);
}

#endif
