/*
 * kernel_axmy_parital_8_dp.cl
 *
 *  Created on: Dec 21, 2011
 *      Author: Matthew Wezowicz
 */

#define BLOCK_SIZE 16
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void matrixAxmyKernelModular8DP(__global double* D, __global double* A, __global double* B,
		__global double* C, int width_A, int width_B, double mod){
	//Geet Workgroup ID
	int bx = get_group_id(0);
	int by = get_group_id(1);

	//Get Local ID
	int tx = get_local_id(0);
	int ty = get_local_id(1);

	//Range of indexies for submatrix of A
	int aBegin= width_A * BLOCK_SIZE * by;
	int aEnd = aBegin + width_A - 1;
	int aStep = BLOCK_SIZE;

	//Range of indecies for sub-matrix of B
	int bBegin = BLOCK_SIZE * bx;
	int bStep = BLOCK_SIZE * width_B;

	//Local storage of sub-matrices of A and B;
	__local double As[BLOCK_SIZE][BLOCK_SIZE];
	__local double Bs[BLOCK_SIZE][BLOCK_SIZE];

	//Temporary storage for result
	double Dsub = 0;

	//Loop over all the sub-maticies of A and B required to compute
	//the result sub-matrix
	for(int a = aBegin, b = bBegin; a < aEnd; a += aStep, b += bStep){
		//Load the matricies from global memory to local memory
		//Each thread loads one element of each sub-matrix
		As[ty][tx] = A[a + width_A * ty + tx];
		Bs[ty][tx] = B[b + width_B * ty + tx];

		//Synchronize threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//Multiply the two sub-matrices together
		for(int i = 0; i < BLOCK_SIZE / 2; i++){
			Dsub += As[ty][i] * Bs[i][tx];
		}
		Dsub = fmod(Dsub, mod);
		for(int i = BLOCK_SIZE / 2; i < BLOCK_SIZE; i++){
			Dsub += As[ty][i] * Bs[i][tx];
		}
		Dsub = fmod(Dsub, mod);

		//Synchronize threads
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	//Calculates the offset inthe result matrix
	int d = width_B * BLOCK_SIZE * by + BLOCK_SIZE * bx;

	//Load, add, and normalize with element from C
	double c = C[d + ty * width_B + tx];
	Dsub = Dsub - c;
	Dsub = fmod((mod + Dsub), mod);

	//Add the sum to the appropriate spot
	D[d + ty * width_B + tx] = Dsub;
}