#pragma once

#ifndef LONG_KERNELS

__kernel void MultipleFinalization_Combine(
	global		const	unsigned int*	lengths,
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localBuffer)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	global float* posterior = posteriorMatrices + tasks[taskId].offset_Pij;

	// calculate combined posterior matrix
	Finalization_CombineMatrices(
		auxiliaryMatrices + tasks[taskId].offset_Aij, 
		posteriorMatrices + tasks[taskId].offset_Pij, 
		lengths[task_i],
		lengths[task_j],
		tasks[taskId].pid,
		localBuffer, 
		posterior, 
		(global float*)0); 
		
	// calculate distance on the basis of score
	if (get_local_id(0) == 0) {
		posterior[0] = 1.0f - localBuffer[0] / min(lengths[task_i], lengths[task_j]);	
		tasks[taskId].numCells = localBuffer[1]; 
	}
}

#else

__kernel void MultipleFinalization_Combine_long(
	global		const	unsigned int*	lengths,
	global				float*			auxiliaryMatrices,
	global				float*			posteriorMatrices,
	global				struct			PosteriorTask* tasks,
						struct			PosteriorSchedule scheduling,
	local				float*			localBuffer,
	global				float*			columns,
				const	int				columnsOffset)
{
	// execute computations
	int taskId = get_group_id(0);
	int task_i = tasks[taskId].i;
	int task_j = tasks[taskId].j;
	
	global float* posterior = posteriorMatrices + tasks[taskId].offset_Pij;

	// calculate combined posterior matrix
	Finalization_CombineMatrices_long(
		auxiliaryMatrices + tasks[taskId].offset_Aij, 
		posteriorMatrices + tasks[taskId].offset_Pij, 
		lengths[task_i],
		lengths[task_j],
		tasks[taskId].pid,
		localBuffer, 
		posterior, 
		(global float*)0,
		columns + taskId * columnsOffset); 
		
	// calculate distance on the basis of score
	if (get_local_id(0) == 0) {
		posterior[0] = 1.0f - localBuffer[0] / min(lengths[task_i], lengths[task_j]);	
		tasks[taskId].numCells = localBuffer[1]; 
	}
}

#endif