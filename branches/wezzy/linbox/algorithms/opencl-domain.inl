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

//#include "helper_functions.cpp"

#include "CL/cl.hpp"

namespace LinBox
{

	/*
	 * ******************************************************
	 * *** Specializations for BlasMatrix<Field> where    ***
	 * *** the Field is Modular<float> or Modular<double> ***
	 * ******************************************************
	 */

	/*
	 * Specialization of addition over
	 * a Modular<double> Field
	 * C = A+B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::add<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !dpKernelsAvailable[0]){
			return BlasMatrixDomainAdd<
				Modular<double>,BlasMatrix<double>,BlasMatrix<double>,BlasMatrix<double> >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.coldim());
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( A.coldim() == C.coldim());
		linbox_check( A.rowdim() == C.rowdim());

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel = dpKernels[0];

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = widthA * heightA;
		globalWorkSize[1] = 1;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
			localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block unitl kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		C = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return C;
	}

	/*
	 * Specialization of addition over
	 * a Modular<float> Field
	 * C = A+B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::add<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !spKernelsAvailable[0]){
			return BlasMatrixDomainAdd<
				Modular<float>,BlasMatrix<float>,BlasMatrix<float>,BlasMatrix<float> >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.coldim());
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( A.coldim() == C.coldim());
		linbox_check( A.rowdim() == C.rowdim());

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel = spKernels[0];

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = widthA * heightA;
		globalWorkSize[1] = 1;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
			localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block unitl kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		C = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return C;
	}

	/*
	 * Specialization of in place addition over
	 * a Modular<double> Field
	 * C += B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::addin<
		BlasMatrix<double>, BlasMatrix<double> >(BlasMatrix<double>& C,
		const BlasMatrix<double>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(C);
		//T = C;
		return add<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,C,B);
	}

	/*
	 * Specialization of in place addition over
	 * a Modular<float> Field
	 * C += B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::addin<
		BlasMatrix<float>, BlasMatrix<float> >(BlasMatrix<float>& C,
		const BlasMatrix<float>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(C);
		//T = C;
		return add<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,C,B);
	}

	/*
	 * Specialization of subtraction over
	 * a Modular<double> Field
	 * C = A-B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::sub<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !dpKernelsAvailable[1]){
			return BlasMatrixDomainSub<
				Modular<double>,BlasMatrix<double>,BlasMatrix<double>,BlasMatrix<double> >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.coldim());
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( A.coldim() == C.coldim());
		linbox_check( A.rowdim() == C.rowdim());

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel = dpKernels[1];

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = widthA * heightA;
		globalWorkSize[1] = 1;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
			localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block unitl kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		C = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return C;
	}

	/*
	 * Specialization of subtraction over
	 * a Modular<float> Field
	 * C = A-B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::sub<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B);

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !spKernelsAvailable[1]){
			return BlasMatrixDomainSub<
				Modular<float>,BlasMatrix<float>,BlasMatrix<float>,BlasMatrix<float> >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.coldim());
		linbox_check( A.rowdim() == B.rowdim());
		linbox_check( A.coldim() == C.coldim());
		linbox_check( A.rowdim() == C.rowdim());

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel = spKernels[1];

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 256;
		localWorkSize[1] = 1;
		globalWorkSize[0] = widthA * heightA;
		globalWorkSize[1] = 1;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
			localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block unitl kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		C = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return C;
	}

	/*
	 * Specialization of in place subtraction over
	 * a Modular<double> Field
	 * C -= B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::subin<
		BlasMatrix<double>, BlasMatrix<double> >(BlasMatrix<double>& C,
		const BlasMatrix<double>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(C);
		//T = C;
		return sub<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,C,B);
	}

	/*
	 * Specialization of in place subtraction over
	 * a Modular<float> Field
	 * C -= B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::subin<
		BlasMatrix<float>, BlasMatrix<float> >(BlasMatrix<float>& C,
		const BlasMatrix<float>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(C);
		//T = C;
		return sub<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,C,B);
	}

	/*
	 * Specialization of Mul for
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
			
		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[2];
		kernelsAvailable &= dpKernelsAvailable[3];
		kernelsAvailable &= dpKernelsAvailable[4];
		kernelsAvailable &= dpKernelsAvailable[5];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMul<
				Modular<double>,BlasMatrix<double>,BlasMatrix<double>,BlasMatrix<double> >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);

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

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		C = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return C;
	}

	/*
	 * Specialization of Mul for
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
		
		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[2];
		kernelsAvailable &= spKernelsAvailable[3];
		kernelsAvailable &= spKernelsAvailable[4];
		kernelsAvailable &= spKernelsAvailable[5];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMul<
				Modular<float>,BlasMatrix<float>,BlasMatrix<float>,BlasMatrix<float> >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);

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

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		C = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return C;
	}

	/*
	 * Specialization of mulin_left for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(A);
		//T = A;
		return mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(A,A,B);
	}

	/*
	 * Specialization of mulin_left for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(A);
		//T = A;
		return mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(A,A,B);
	}

	/*
	 * Specialization of mulin_right for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(B);
		//T = B;
		return mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(B,A,B);
	}

	/*
	 * Specialization of mulin_right for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(B);
		//T = B;
		return mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(B,A,B);
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<double> Field
	 * D = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::muladd<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& D, const double& beta, const BlasMatrix<double>& C,
		const double& alpha, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[6];
		kernelsAvailable &= dpKernelsAvailable[7];
		kernelsAvailable &= dpKernelsAvailable[8];
		kernelsAvailable &= dpKernelsAvailable[9];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<double>, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >()(
				_F,D,beta,C,alpha,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		double p = _F.characteristic();
		double tempAlpha = fmod(alpha, p);
		double tempBeta = fmod(beta, p);

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[9];
		}
		else if(p <=(1<<24)){
			selectedKernel = dpKernels[8];
		}
		else if(p <=(1<<25)){
			selectedKernel = dpKernels[7];
		}
		else{
			selectedKernel = dpKernels[6];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_double), (void*)&tempAlpha);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_double), (void*)&tempBeta);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 7, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 8, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //does not work because of const being used pointlessly

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<float> Field
	 * D = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::muladd<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& D, const float& beta, const BlasMatrix<float>& C,
		const float& alpha, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[6];
		kernelsAvailable &= spKernelsAvailable[7];
		kernelsAvailable &= spKernelsAvailable[8];
		kernelsAvailable &= spKernelsAvailable[9];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<float>, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >()(
				_F,D,beta,C,alpha,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		float p = _F.characteristic();
		float tempAlpha = fmod(alpha, p);
		float tempBeta = fmod(beta, p);

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[9];
		}
		else if(p <=(1<<9)){
			selectedKernel = spKernels[8];
		}
		else if(p <=(1<<10)){
			selectedKernel = spKernels[7];
		}
		else{
			selectedKernel = spKernels[6];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_float), (void*)&tempAlpha);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_float), (void*)&tempBeta);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 7, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 8, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //does not work because of const being used pointlessly

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<double> Field
	 * Places the results into the first genreral dense matrix
	 * C = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::muladdin<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		const double& beta, BlasMatrix<double>& C, const double& alpha,
		const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(C);
		//T = C;
		return muladd<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,beta,C,alpha,A,B);
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<float> Field
	 * Places the results into the first genreral dense matrix
	 * C = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::muladdin<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		const float& beta, BlasMatrix<float>& C, const float& alpha,
		const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(C);
		//T = C;
		return muladd<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,beta,C,alpha,A,B);
	}

	/*
	 * Specialization of multiplication with scaling over
	 * a Modular<double> Field
	 * C = alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::mul<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const double& alpha, const BlasMatrix<double>& A,
		const BlasMatrix<double>& B) const{

		return muladd<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,0,C,alpha,A,B);
	}

	/*
	 * Specialization of multiplication with scaling over
	 * a Modular<float> Field
	 * C = alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::mul<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const float& alpha, const BlasMatrix<float>& A,
		const BlasMatrix<float>& B) const{

		return muladd<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,0,C,alpha,A,B);
	}

	/*
	 * Specialization of apxy for
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
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[10];
		kernelsAvailable &= dpKernelsAvailable[11];
		kernelsAvailable &= dpKernelsAvailable[12];
		kernelsAvailable &= dpKernelsAvailable[13];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !kernelsAvailable){
			D = mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B);
			return addin<BlasMatrix<double>, BlasMatrix<double> >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[13];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[12];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[11];
		}
		else{
			selectedKernel = dpKernels[10];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of apxy for
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
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[10];
		kernelsAvailable &= spKernelsAvailable[11];
		kernelsAvailable &= spKernelsAvailable[12];
		kernelsAvailable &= spKernelsAvailable[13];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !kernelsAvailable){
			D = mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B);
			return addin<BlasMatrix<float>, BlasMatrix<float> >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[13];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[12];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[11];
		}
		else{
			selectedKernel = spKernels[10];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of apxyin for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(C);
		//T = C;
		return axpy<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B,C);
	}

	/*
	 * Specialization of apxyin for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(C);
		//T = C;
		return axpy<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B,C);
	}

	/*
	 * Specialization of maxpy for
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
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[14];
		kernelsAvailable &= dpKernelsAvailable[15];
		kernelsAvailable &= dpKernelsAvailable[16];
		kernelsAvailable &= dpKernelsAvailable[17];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !kernelsAvailable){
			BlasMatrix<double> T(D.rowdim(),D.coldim());
			T = mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(T,A,B);
			return sub<BlasMatrix<double>, BlasMatrix<double> >(D,C,T);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[17];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[16];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[15];
		}
		else{
			selectedKernel = dpKernels[14];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of maxpy for
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
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[14];
		kernelsAvailable &= spKernelsAvailable[15];
		kernelsAvailable &= spKernelsAvailable[16];
		kernelsAvailable &= spKernelsAvailable[17];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !kernelsAvailable){
			BlasMatrix<float> T(D.rowdim(),D.coldim());
			T = mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(T,A,B);
			return sub<BlasMatrix<float>, BlasMatrix<float> >(D,C,T);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[17];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[16];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[15];
		}
		else{
			selectedKernel = spKernels[14];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of maxpyin for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(C);
		//T = C;
		return maxpy<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B,C);
	}

	/*
	 * Specialization of maxpyin for
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

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(C);
		//T = C;
		return maxpy<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B,C);
	}

	/*
	 * Specialization of axmy for
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
			cl_double, BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[18];
		kernelsAvailable &= dpKernelsAvailable[19];
		kernelsAvailable &= dpKernelsAvailable[20];
		kernelsAvailable &= dpKernelsAvailable[21];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !doubleSupported || !kernelsAvailable){
			D = mul<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(D,A,B);
			return subin<BlasMatrix<double>, BlasMatrix<double> >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<double> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<double> >(C);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[21];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[20];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[19];
		}
		else{
			selectedKernel = dpKernels[18];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_double), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<double> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of axmy for
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
			cl_float, BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B,C);
			
		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[18];
		kernelsAvailable &= spKernelsAvailable[19];
		kernelsAvailable &= spKernelsAvailable[20];
		kernelsAvailable &= spKernelsAvailable[21];

		//If it is not capable or not setup properly use default implementation
		if(!memLevelsAllowed || !setupCorrect || !kernelsAvailable){
			D = mul<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(D,A,B);
			return subin<BlasMatrix<float>, BlasMatrix<float> >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<float> >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<float> >(C);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[21];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[20];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[19];
		}
		else{
			selectedKernel = spKernels[18];
		}

		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int widthA = ((A.coldim() + 15) / 16) * 16;
		int heightA = ((A.rowdim() + 15) / 16) * 16;
		int widthB = ((B.coldim() + 15) / 16) * 16;

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_float), (void*)&p);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Read back buffer
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<float> >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return D;
	}

	/*
	 * Specialization of axmyin for
	 * multiplying two general dense matrices
	 * and subtracts a third general dense matrix from it
	 * over a Modular<double> Field.
	 * Places the results into the first general dense matrix
	 * C = A*B - C
	 */
	template <>
	template <>
	BlasMatrix<double>& OpenCLMatrixDomain<Modular<double> >::axmyin<
		BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(
		BlasMatrix<double>& C, const BlasMatrix<double>& A, const BlasMatrix<double>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<double> T(C);
		//T = C;
		return axmy<BlasMatrix<double>, BlasMatrix<double>, BlasMatrix<double> >(C,A,B,C);
	}

	/*
	 * Specialization of axmyin for
	 * multiplying two general dense matrices
	 * and subtracts a third general dense matrix from it
	 * over a Modular<float> Field.
	 * Places the results into the first general dense matrix
	 * C = A*B - C
	 */
	template <>
	template <>
	BlasMatrix<float>& OpenCLMatrixDomain<Modular<float> >::axmyin<
		BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(
		BlasMatrix<float>& C, const BlasMatrix<float>& A, const BlasMatrix<float>& B) const{

		//Create a temporary matrix for psuedo inline
		//BlasMatrix<float> T(C);
		//T = C;
		return axmy<BlasMatrix<float>, BlasMatrix<float>, BlasMatrix<float> >(C,A,B,C);
	}

} //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_INL

