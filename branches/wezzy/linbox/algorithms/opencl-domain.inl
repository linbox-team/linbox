/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain.inl
 * Copyright (C) 2011 Matthew Wezowicz, David Saunders
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_opencl_matrix_domain_INL
#define __LINBOX_opencl_matrix_domain_INL

#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/matrix/factorized-matrix.h"

#include "CL/cl.hpp"
//#include "helper_functions.hpp" -- For debugging only

namespace LinBox
{

	/*
	 * ******************************************************
	 * *** Specializations for BlasMatrix<Field> where    ***
	 * *** the Field is Modular<float> or Modular<double> ***
	 * ******************************************************
	 */

	/* Specialization of Mul for
	 * multiplying two general dense matrices
	 * over a Modular<double> Field.
	 * C = A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::mul<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect && !doubleSupported){
			return BlasMatrixDomainMul<
				Modular<double>,BlasMatrix<double>,BlasMatrix<double>,BlasMatrix<double> >()(_F,C,A,B);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());

		cl_mem bufferC = createMatrixBuffer<cl_double, BlasMatrix<double> >(C);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[5];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[4];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[3];
		}
		else{
			selectedKernel = dpKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		C = readMatrixBuffer<cl_double, BlasMatrix<double> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return C;
	}

	/* Specialization of Mul for
	 * multiplying two general dense matrices
	 * over a Modular<float> Field.
	 * C = A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::mul<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect){
			return BlasMatrixDomainMul<
				Modular<float>,BlasMatrix<float>,BlasMatrix<float>,BlasMatrix<float> >()(_F,C,A,B);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());

		cl_mem bufferC = createMatrixBuffer<cl_float, BlasMatrix<float> >(C);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[5];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[4];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[3];
		}
		else{
			selectedKernel = spKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		C = readMatrixBuffer<cl_float, BlasMatrix<float> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return C;
	}

	/* Specialization of mulin_left for
	 * multiplying two general dense matrices
	 * over a Modular<double> Field.
	 * Places result into the left matrix.
	 * A = A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::mulin_left<
		BlasMatrix<double>, BlasMatrix<double> >( BlasMatrix<double>& A,
		const BlasMatrix<double>& B) const{

		BlasMatrix<double> T(A);
		T = A;
		return mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(A,T,B);
	}

	/* Specialization of mulin_left for
	 * multiplying two general dense matrices
	 * over a Modular<float> Field.
	 * Places the result into the left matrix.
	 * A = A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::mulin_left<
		BlasMatrix<float>, BlasMatrix<float> >( BlasMatrix<float>& A,
		const BlasMatrix<float>& B) const{

		BlasMatrix<float> T(A);
		T = A;
		return mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(A,T,B);
	}

	/* Specialization of mulin_right for
	 * multiplying two general dense matrices
	 * over a Modular<double> Field.
	 * Places the result into the right matrix.
	 * B = A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::mulin_right<
		BlasMatrix<double>, BlasMatrix<double> >(const BlasMatrix<double>& A,
		BlasMatrix<double>& B) const{

		BlasMatrix<double> T(B);
		T = B;
		return mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(B,A,T);
	}

	/* Specialization of mulin_right for
	 * multiplying two general dense matrices
	 * over a Modular<float> Field.
	 * Places the result into the right matrix.
	 * B = A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::mulin_right<
		BlasMatrix<float>, BlasMatrix<float> >(const BlasMatrix<float>& A,
		BlasMatrix<float>& B) const{

		BlasMatrix<float> T(B);
		T = B;
		return mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(B,A,T);
	}

	/* Specialization of apxy for
	 * multiplying two general dense matrices
	 * and adding a third general dense matrix
	 * over a Modular<double> Field.
	 * D = A*B + C
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::axpy<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& D, const BlasMatrix<double>& A, const BlasMatrix<double>& B,
		const BlasMatrix<double>& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C,D);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect && !doubleSupported){
			D = mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B);
			return addin<BlasMatrix<double>, BlasMatrix<double> >(D,C);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		cl_mem bufferT = createMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[5];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[4];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[3];
		}
		else{
			selectedKernel = dpKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass 1st(mul) kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch 1st(mul) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		cl_mem bufferD = createMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferC = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		//Select modular addition kernel
		selectedKernel = dpKernels[0];

		//Pass 2nd(add) kernel arguments
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = (heightA * widthB);
		globalWorkSize[1] = 1;

		//Launch 2nd(add) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		D = readMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferT);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return D;
	}

	/* Specialization of apxy for
	 * multiplying two general dense matrices
	 * and adding a third general dense matrix
	 * over a Modular<float> Field.
	 * D = A*B + C
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::axpy<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& D, const BlasMatrix<float>& A, const BlasMatrix<float>& B,
		const BlasMatrix<float>& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C,D);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect){
			D = mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B);
			return addin<BlasMatrix<float>, BlasMatrix<float> >(D,C);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		cl_mem bufferT = createMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[5];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[4];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[3];
		}
		else{
			selectedKernel = spKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass 1st(mul) kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch 1st(mul) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		cl_mem bufferD = createMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferC = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		//Select modular addition kernel
		selectedKernel = dpKernels[0];

		//Pass 2nd(add) kernel arguments
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = (heightA * widthB);
		globalWorkSize[1] = 1;

		//Launch 2nd(add) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		D = readMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferT);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return D;
	}

	/* Specialization of apxyin for
	 * multiplying two general dense matrices
	 * and adding a third general dense matrix
	 * over a Modular<double> Field.
	 * Places the result into the first matrix.
	 * C += A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::axpyin<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		BlasMatrix<double> T(C);
		T = C;
		return axpy<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B,T);
	}

	/* Specialization of apxyin for
	 * multiplying two general dense matrices
	 * and adding a third general dense matrix
	 * over a Modular<float> Field.
	 * Places the result into the first matrix.
	 * C += A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::axpyin<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		BlasMatrix<float> T(C);
		T = C;
		return axpy<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B,T);
	}

	/* Specialization of maxpy for
	 * multiplying two general dense matrices
	 * and subtracts it from a third general dense matrix
	 * over a Modular<double> Field.
	 * D = C - A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::maxpy<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& D, const BlasMatrix<double>& A, const BlasMatrix<double>& B,
		const BlasMatrix<double>& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C,D);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect && !doubleSupported){
			BlasMatrix<double> T(D.rowdim(),D.coldim());
			T = mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(T,A,B);
			return sub<BlasMatrix<double>, BlasMatrix<double> >(D,C,T);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		cl_mem bufferT = createMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[5];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[4];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[3];
		}
		else{
			selectedKernel = dpKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass 1st(mul) kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch 1st(mul) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		cl_mem bufferD = createMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferC = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		//Select modular subtraction kernel
		selectedKernel = dpKernels[1];

		//Pass 2nd(sub) kernel arguments
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = (heightA * widthB);
		globalWorkSize[1] = 1;

		//Launch 2nd(sub) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		D = readMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferT);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return D;
	}

	/* Specialization of maxpy for
	 * multiplying two general dense matrices
	 * and subtracts it from a third general dense matrix
	 * over a Modular<float> Field.
	 * D = C - A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::maxpy<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& D, const BlasMatrix<float>& A, const BlasMatrix<float>& B,
		const BlasMatrix<float>& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C,D);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect){
			BlasMatrix<float> T(D.rowdim(),D.coldim());
			T = mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(T,A,B);
			return sub<BlasMatrix<float>, BlasMatrix<float> >(D,C,T);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		cl_mem bufferT = createMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[5];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[4];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[3];
		}
		else{
			selectedKernel = spKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass 1st(mul) kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch 1st(mul) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		cl_mem bufferD = createMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferC = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		//Select modular subtraction kernel
		selectedKernel = dpKernels[1];

		//Pass 2nd(sub) kernel arguments
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = (heightA * widthB);
		globalWorkSize[1] = 1;

		//Launch 2nd(sub) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		D = readMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferT);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return D;
	}

	/* Specialization of maxpyin for
	 * multiplying two general dense matrices
	 * and subtracts it from a third general dense matrix
	 * over a Modular<double> Field.
	 * Places the results into the first gernal dense matrix.
	 * C -= A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::maxpyin<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		BlasMatrix<double> T(C);
		T = C;
		return maxpy<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B,T);
	}

	/* Specialization of maxpyin for
	 * multiplying two general dense matrices
	 * and subtracts it from a third general dense matrix
	 * over a Modular<float> Field.
	 * Places the results into the first gernal dense matrix.
	 * C -= A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::maxpyin<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		BlasMatrix<float> T(C);
		T = C;
		return maxpy<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B,T);
	}

	/* Specialization of axmy for
	 * multiplying two general dense matrices
	 * and subtracts a third general dense matrix from it
	 * over a Modular<double> Field.
	 * D = A*B - C
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::axmy<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& D, const BlasMatrix<double>& A, const BlasMatrix<double>& B,
		const BlasMatrix<double>& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C,D);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect && !doubleSupported){
			D = mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B);
			return subin<BlasMatrix<double>, BlasMatrix<double> >(D,C);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		cl_mem bufferT = createMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[5];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[4];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[3];
		}
		else{
			selectedKernel = dpKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass 1st(mul) kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch 1st(mul) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		cl_mem bufferD = createMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferC = createAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		//Select modular subtraction kernel
		selectedKernel = dpKernels[1];

		//Pass 2nd(sub) kernel arguments
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = (heightA * widthB);
		globalWorkSize[1] = 1;

		//Launch 2nd(sub) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		D = readMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferT);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return D;
	}

	/* Specialization of mapxy for
	 * multiplying two general dense matrices
	 * and subtracts a third general dense matrix from it
	 * over a Modular<float> Field.
	 * D = A*B - C
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::axmy<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& D, const BlasMatrix<float>& A, const BlasMatrix<float>& B,
		const BlasMatrix<float>& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C,D);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed && !setupCorrect){
			D = mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B);
			return subin<BlasMatrix<float>, BlasMatrix<float> >(D,C);
		}

		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		cl_mem bufferT = createMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[5];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[4];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[3];
		}
		else{
			selectedKernel = spKernels[2];
		}

		int widthA = ((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16;
		int heightA = ((A.rowdim() / 16) + (A.rowdim() % 16 == 0 ? 0 : 1)) * 16;
		int widthB = ((B.coldim() / 16) + (B.coldim() % 16 == 0 ? 0 : 1)) * 16;

		//Pass 1st(mul) kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch 1st(mul) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		cl_mem bufferD = createMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferC = createAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		//Select modular subtraction kernel
		selectedKernel = dpKernels[1];

		//Pass 2nd(sub) kernel arguments
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferT);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Set NDRange
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = (heightA * widthB);
		globalWorkSize[1] = 1;

		//Launch 2nd(sub) kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
				localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		D = readMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferT);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const being used pointlessly

		return D;
	}

	/* Specialization of axmyin for
	 * multiplying two general dense matrices
	 * and subtracts a third general dense matrix from it
	 * over a Modular<double> Field.
	 * PLaces the results into the first general dense matrix
	 * C = A*B - C
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::axmyin<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		BlasMatrix<double> T(C);
		T = C;
		return axmy<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B,T);
	}

	/* Specialization of axmyin for
	 * multiplying two general dense matrices
	 * and subtracts a third general dense matrix from it
	 * over a Modular<float> Field.
	 * PLaces the results into the first general dense matrix
	 * C = A*B - C
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::axmyin<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		BlasMatrix<float> T(C);
		T = C;
		return axmy<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B,T);
	}

} //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_INL

