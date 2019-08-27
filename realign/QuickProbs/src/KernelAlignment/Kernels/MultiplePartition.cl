#pragma once

#ifndef LONG_KERNELS

__kernel void MultiplePartition_ComputeForward(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			PartitionFunctionParams* globalParams MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				partition_t*	localMemory)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PARTITION_PARAMS
		// copy scoring matrices from global to local memory
		local struct PartitionFunctionParams localParams;
		local struct PartitionFunctionParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct PartitionFunctionParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Partition_ComputeForward(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij);
}


__kernel void MultiplePartition_ComputeReverse(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			PartitionFunctionParams* globalParams MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				partition_t*	localMemory)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PARTITION_PARAMS
		// copy scoring matrices from global to local memory
		local struct PartitionFunctionParams localParams;
		local struct PartitionFunctionParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct PartitionFunctionParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Partition_ComputeReverse(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij,
		posteriorMatrices + tasks[taskId].offset_Pij);
}



#else


__kernel void MultiplePartition_ComputeForward_long(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			PartitionFunctionParams* globalParams MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				partition_t*	localMemory,
	global				float*			columns,
				const	int				columnsOffset)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PARTITION_PARAMS
		// copy scoring matrices from global to local memory
		local struct PartitionFunctionParams localParams;
		local struct PartitionFunctionParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct PartitionFunctionParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Partition_ComputeForward_long(
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


__kernel void MultiplePartition_ComputeReverse_long(
	global		const	char*			sequences,
	global		const	unsigned int*	offsets,
	global		const	unsigned int*	lengths,
	global		const	struct			PartitionFunctionParams* globalParams MAX_SIZE(sizeof(struct PartitionFunctionParams)),
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				partition_t*	localMemory,
	global				float*			columns,
				const	int				columnsOffset)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	#ifdef COPY_PARTITION_PARAMS
		// copy scoring matrices from global to local memory
		local struct PartitionFunctionParams localParams;
		local struct PartitionFunctionParams *params = &localParams;

		event_t eventCopied = async_work_group_copy(
				(local int*)params, (global int*)globalParams, sizeof(struct PartitionFunctionParams) / sizeof(int), 0);
		wait_group_events(1, &eventCopied);
	#else
		#define params globalParams
	#endif

	Partition_ComputeReverse_long(
		sequences + offsets[task_i],
		sequences + offsets[task_j],
		lengths[task_i],
		lengths[task_j],
		scheduling,
		params,
		localMemory,
		auxiliaryMatrices + tasks[taskId].offset_Aij,
		posteriorMatrices + tasks[taskId].offset_Pij,
		columns + taskId * columnsOffset);
}

#endif