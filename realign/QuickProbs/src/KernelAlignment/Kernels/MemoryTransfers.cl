// This file contains some memory transfers definitions
#pragma once

#define swap(x,y,temp) {temp = (x); x = (y); y = (temp);}

/// <summary>
/// Copies num elements starting from offset from src to dst location.
/// </summary>
#define __memcpy(dst,src,num,offset) { for (int i = offset; i < offset + num; i++) { dst[i] = src[i]; } }

/// <summary>
/// Sets num elements to value starting from ptr location.
/// </summary>
#define __memset(ptr,value,num) { for (int i = 0; i < num; i++) { (ptr)[i] = (value); } }

void memcpy(local float* dst, global float* src, int num)
{
	__memcpy(dst, src, num, 0);
}

void group_multiply(global float* ptr, int size, float value)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0)) { ptr[i] *= value; }
}

void group_copy_g2g(global float* dst, global float* src, int size)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0)) { dst[i] = src[i]; }
}

void group_copy_g2l(local float* dst, global float* src, int size)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0)) { dst[i] = src[i]; }
}

void group_copy_l2g(global float* dst, local float* src, int size)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0)) { dst[i] = src[i]; }
}

void group_weightcopy_g2l(local float* dst, global float* src, float weight, int size)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0)) { dst[i] = src[i] * weight; }
}

void group_weightcopy_l2g(global float* dst, local float* src, float weight, int size)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0)) { dst[i] = src[i] * weight; }
}

void group_memset(global float* dst, int size, float value)
{
	for (int i = get_local_id(0); i < size; i += get_local_size(0))
	{
		dst[i] = value;
	}
} 


