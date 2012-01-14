/*
 * kernel_partial_1024_dp.cl
 *
 *  Created on: Jul 13, 2011
 *      Author: Matthew Wezowicz
 */

#define BLOCK_SIZE 16
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void matrixMulKernelModular1024DP(__global double* C, __global double* A, __global double* B,
		const int widthA, const int widthB, const double mod){
	//Get Workgroup ID
	int bx = get_group_id(0);
	int by = get_group_id(1);

	//Get Local ID
	int tx = get_local_id(0);
	int ty = get_local_id(1);

	//Range of indecies for sub-matrix of A
	int aBegin = widthA * BLOCK_SIZE * by;
	int aEnd = aBegin + widthA - 1;
	int aStep = BLOCK_SIZE;

	//Range of indecies for sub-matrix of B
	int bBegin = BLOCK_SIZE * bx;
	int bStep = BLOCK_SIZE * widthB;

	//Local storage of sub-matrices of A and B
	__local double As[BLOCK_SIZE][BLOCK_SIZE];
	__local double Bs[BLOCK_SIZE][BLOCK_SIZE];

	//Temporary storage for result
	double Csub = 0;

	//Counter for modulus every 32 iterations
	int mCount = 0;

	//Loop over all the sub-matrices of A and B required to compute
	//the result sub-matrix
	for(int a = aBegin, b = bBegin; a < aEnd; a += aStep, b += bStep){
		//Load the matrices from global memory to local memory
		//Each thread loads one element of each sub-matrix
		As[ty][tx] = A[a + widthA * ty + tx];
		Bs[ty][tx] = B[b + widthB * ty + tx];

		//Synchronize threads
		barrier(CLK_LOCAL_MEM_FENCE);

		//Multiply the two sub-matrices together
		Csub += As[ty][0] * Bs[0][tx];
		Csub += As[ty][1] * Bs[1][tx];
		Csub += As[ty][2] * Bs[2][tx];
		Csub += As[ty][3] * Bs[3][tx];
		Csub += As[ty][4] * Bs[4][tx];
		Csub += As[ty][5] * Bs[5][tx];
		Csub += As[ty][6] * Bs[6][tx];
		Csub += As[ty][7] * Bs[7][tx];
		Csub += As[ty][8] * Bs[8][tx];
		Csub += As[ty][9] * Bs[9][tx];
		Csub += As[ty][10] * Bs[10][tx];
		Csub += As[ty][11] * Bs[11][tx];
		Csub += As[ty][12] * Bs[12][tx];
		Csub += As[ty][13] * Bs[13][tx];
		Csub += As[ty][14] * Bs[14][tx];
		Csub += As[ty][15] * Bs[15][tx];

		mCount++;

		//fmod every 1024 iterations
		if(mCount == 64){
			Csub = fmod(Csub, mod);
			mCount = 0;
		}

		//Synchronize threads
		barrier(CLK_LOCAL_MEM_FENCE);
	}
	//Calls fmod once to normalize the sum
	Csub = fmod(Csub, mod);

	//Calculates the offset in the result matrix and add the sum to the
	//appropriate spot
	int c = widthB * BLOCK_SIZE * by + BLOCK_SIZE * bx;
	C[c + ty * widthB + tx] = Csub;
}
