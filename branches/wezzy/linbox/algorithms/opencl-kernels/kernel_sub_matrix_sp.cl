/*
 * kernel_sub_matrix_sp.cl
 *
 *  Created on: Jul 18, 2011
 *      Author: Matthew Wezowicz
 */

#define BLOCK_SIZE 256

__kernel void vector_sum_kernel(__global float* C, __global float* A, __global float* B, float mod){
	//Get Workgroup ID
	int bx = get_group_id(0);

	//Get Local ID
	int tx = get_local_id(0);

	//Calculate global ID
	int gx = bx * BLOCK_SIZE + tx;

	//Load from matrices
	float a = A[gx];
	float b = B[gx];

	//Sum and fmod()
	float c = fmod((a - b), mod);

	//Write to result matrix
	C[gx] = c;
}
