/*
 * kernel_add_matrix_dp.cl
 *
 *  Created on: Jul 18, 2011
 *      Author: Matthew Wezowicz
 */

#define BLOCK_SIZE 256
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void vector_sum_kernel(__global double* C, __global double* A, __global double* B, double mod){
	//Get Workgroup ID
	int bx = get_group_id(0);

	//Get Local ID
	int tx = get_local_id(0);

	//Calculate global ID
	int gx = bx * BLOCK_SIZE + tx;

	//Load from matrices
	double a = A[gx];
	double b = B[gx];

	//Sum and fmod()
	double c = fmod((a + b), mod);

	//Write to result matrix
	C[gx] = c;
}
