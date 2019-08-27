#pragma once

#ifndef LONG_KERNELS

/// <summary>
/// </summary>
/// <param name=lengths></param>
/// <param name=tasks></param>
/// <param name=scheduling></param>
/// <param name=sparseMatrices></param>
/// <param name=posteriorMatrices></param>
/// <returns></returns>
__kernel void MultipleSparse_Compute(
	global	const	unsigned int* lengths,
	global			float* sparseMatrices,
	global			float* posteriorMatrices,
	global	const	struct PosteriorTask* tasks,
					struct PosteriorSchedule scheduling,
	local			float* localMemory,
	int				transpose)
{
	int taskId = get_group_id(0);
	int i = tasks[taskId].i;
	int j = tasks[taskId].j;

	if (!transpose) {
		SparseMatrix_Generate(
			lengths[i],
			lengths[j],
			posteriorMatrices + tasks[taskId].offset_Pij,
			localMemory,
			sparseMatrices + tasks[taskId].offset_Aij);
	} else {
		SparseMatrix_GenerateTranspose(
			lengths[i],
			lengths[j],
			posteriorMatrices + tasks[taskId].offset_Pij,
			localMemory,
			sparseMatrices+ tasks[taskId].offset_Aij);
	}
}

#else

/// <summary>
/// </summary>
/// <param name=lengths></param>
/// <param name=tasks></param>
/// <param name=scheduling></param>
/// <param name=sparseMatrices></param>
/// <param name=posteriorMatrices></param>
/// <returns></returns>
__kernel void MultipleSparse_Compute_long(
	global	const	unsigned int* lengths,
	global			float* sparseMatrices,
	global			float* posteriorMatrices,
	global	const	struct PosteriorTask* tasks,
					struct PosteriorSchedule scheduling,
	local			float* localMemory,
	int				transpose)
{
	int taskId = get_group_id(0);
	int i = tasks[taskId].i;
	int j = tasks[taskId].j;

	if (!transpose) {
		SparseMatrix_Generate_verylong(
			lengths[i],
			lengths[j],
			posteriorMatrices + tasks[taskId].offset_Pij,
			localMemory,
			sparseMatrices + tasks[taskId].offset_Aij);
	} else {
		SparseMatrix_GenerateTranspose_long(
			lengths[i],
			lengths[j],
			posteriorMatrices + tasks[taskId].offset_Pij,
			localMemory,
			sparseMatrices+ tasks[taskId].offset_Aij);
	}
}

#endif