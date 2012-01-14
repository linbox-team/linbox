/*
 * kernel_maxpy_parital_32_dp.cl
 *
 *  Created on: Dec 21, 2011
 *      Author: Matthew Wezowicz
 */

#define BLOCK_SIZE 16
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void matrixMaxpyKernelModular32DP(__global double* D, __global double* A, __global double* B,
		__global double* C, const int widthA, const int widthB, const double mod){
	//Geet Workgroup ID
	int bx = get_group_id(0);
	int by = get_group_id(1);

	//Get Local ID
	int tx = get_local_id(0);
	int ty = get_local_id(1);

	//Range of indexies for submatrix of A
	int aBegin= widthA * BLOCK_SIZE * by;
	int aEnd = aBegin + widthA - 1;
	int aStep = BLOCK_SIZE;

	//Range of indecies for sub-matrix of B
	int bBegin = BLOCK_SIZE * bx;
	int bStep = BLOCK_SIZE * widthB;

	//Local storage of sub-matrices of A and B;
	__local double As[BLOCK_SIZE][BLOCK_SIZE];
	__local double Bs[BLOCK_SIZE][BLOCK_SIZE];

	//Temporary storage for result
	double Dsub = 0;

	//Counter for modulus every 32 iterations
	int mCount = 0;

	//Loop over all the sub-maticies of A and B required to compute
	//the result sub-matrix
	for(int a = aBegin, b = bBegin; a < aEnd; a += aStep, b += bStep){
		//Load the matricies from global memory to local memory
		//Each thread loads one element of each sub-matrix
		As[ty][tx] = A[a + widthA * ty + tx];
		Bs[ty][tx] = B[b + widthB * ty + tx];

		//Synchronize threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//Multiply the two sub-matrices together
		Dsub += As[ty][0] * Bs[0][tx];
		Dsub += As[ty][1] * Bs[1][tx];
		Dsub += As[ty][2] * Bs[2][tx];
		Dsub += As[ty][3] * Bs[3][tx];
		Dsub += As[ty][4] * Bs[4][tx];
		Dsub += As[ty][5] * Bs[5][tx];
		Dsub += As[ty][6] * Bs[6][tx];
		Dsub += As[ty][7] * Bs[7][tx];
		Dsub += As[ty][8] * Bs[8][tx];
		Dsub += As[ty][9] * Bs[9][tx];
		Dsub += As[ty][10] * Bs[10][tx];
		Dsub += As[ty][11] * Bs[11][tx];
		Dsub += As[ty][12] * Bs[12][tx];
		Dsub += As[ty][13] * Bs[13][tx];
		Dsub += As[ty][14] * Bs[14][tx];
		Dsub += As[ty][15] * Bs[15][tx];

		mCount++;

		//fmod every 32 iterations
		if(mCount == 2){
			Dsub = fmod(Dsub, mod);
			mCount = 0;
		}

		//Synchronize threads
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	//Calls fmod once to normalize the sum
	Dsub = fmod(Dsub, mod);

	//Calculates the offset inthe result matrix
	int d = widthB * BLOCK_SIZE * by + BLOCK_SIZE * bx;

	//Load, add, and normalize with element from C
	double c = C[d + ty * widthB + tx];
	Dsub = c - Dsub;
	Dsub = fmod((mod + Dsub), mod);

	//Add the sum to the appropriate spot
	D[d + ty * widthB + tx] = Dsub;
}