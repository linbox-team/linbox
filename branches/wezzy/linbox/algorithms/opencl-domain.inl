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

#include <cstdio>
#include <pthread.h>
#include "linbox/matrix/blas-matrix.h"

#include "/home/mwezz/timer.h"

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
	 * Specialization of Mul for
	 * multiplying two general dense matrices
	 * over a Modular<double> Field.
	 * C = A*B
	 */
	template <>
	template <>
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::mul<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& C, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B) const{

		OpenCLTimer timer;
		double duration = 0;
		timer.tic();
		
		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(C,A,B);

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[0];
		kernelsAvailable &= dpKernelsAvailable[1];
		kernelsAvailable &= dpKernelsAvailable[2];
		kernelsAvailable &= dpKernelsAvailable[3];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMul<
				Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(_F,C,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			return BlasMatrixDomainMul<
				Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		
		duration = timer.toc();
		printf("Checks: %lf\n", duration);
		timer.tic();

		//Lock the device
		pthread_mutex_lock(deviceLock);

		timer.tic();
		
		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(B);
		
		duration = timer.toc();
		printf("Matrix Transfer: %lf\n", duration);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[3];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[2];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[1];
		}
		else{
			selectedKernel = dpKernels[0];
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

		timer.tic();
		
		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(commandQue, selectedKernel, 2, NULL, globalWorkSize,
			localWorkSize, 0, NULL, NULL);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Block until kernel finishes
		tempErrcode = clFinish(commandQue);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually
		
		duration = timer.toc();
		printf("Execution: %lf\n", duration);
		timer.tic();

		//Read back buffer
		C = oclReadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(bufferC, C);
		
		duration = timer.toc();
		printf("Read back: %lf\n", duration);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);
		
		duration = timer.toc();
		printf("Overall time: %lf\n", duration);

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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::mul<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& C, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(C,A,B);

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[0];
		kernelsAvailable &= spKernelsAvailable[1];
		kernelsAvailable &= spKernelsAvailable[2];
		kernelsAvailable &= spKernelsAvailable[3];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMul<
				Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(_F,C,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			return BlasMatrixDomainMul<
				Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(_F,C,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferC = oclCreateMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(C);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[3];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[2];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[1];
		}
		else{
			selectedKernel = spKernels[0];
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
		C = oclReadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(bufferC, C);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferC);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::mulin_left<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >( BlasMatrix<Modular<double> >& A,
		const BlasMatrix<Modular<double> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[0];
		kernelsAvailable &= dpKernelsAvailable[1];
		kernelsAvailable &= dpKernelsAvailable[2];
		kernelsAvailable &= dpKernelsAvailable[3];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulin<Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(_F,A,B);
		}

		return mul<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(A,A,B);
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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::mulin_left<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >( BlasMatrix<Modular<float> >& A,
		const BlasMatrix<Modular<float> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[0];
		kernelsAvailable &= spKernelsAvailable[1];
		kernelsAvailable &= spKernelsAvailable[2];
		kernelsAvailable &= spKernelsAvailable[3];


		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulin<Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(_F,A,B);
		}

		return mul<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(A,A,B);
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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::mulin_right<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(const BlasMatrix<Modular<double> >& A,
		BlasMatrix<Modular<double> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[0];
		kernelsAvailable &= dpKernelsAvailable[1];
		kernelsAvailable &= dpKernelsAvailable[2];
		kernelsAvailable &= dpKernelsAvailable[3];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulin<Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(_F,A,B);
		}

		return mul<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(B,A,B);
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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::mulin_right<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(const BlasMatrix<Modular<float> >& A,
		BlasMatrix<Modular<float> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[0];
		kernelsAvailable &= spKernelsAvailable[1];
		kernelsAvailable &= spKernelsAvailable[2];
		kernelsAvailable &= spKernelsAvailable[3];


		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulin<Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(_F,A,B);
		}

		return mul<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(B,A,B);
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<double> Field
	 * D = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::muladd<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& D, const double& beta, const BlasMatrix<Modular<double> >& C,
		const double& alpha, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[4];
		kernelsAvailable &= dpKernelsAvailable[5];
		kernelsAvailable &= dpKernelsAvailable[6];
		kernelsAvailable &= dpKernelsAvailable[7];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<double>, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >()(
				_F,D,beta,C,alpha,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			return BlasMatrixDomainMulAdd<
				Modular<double>, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >()(
				_F,D,beta,C,alpha,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(C);

		double p = _F.characteristic();
		double tempAlpha = fmod(alpha, p);
		double tempBeta = fmod(beta, p);

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[7];
		}
		else if(p <=(1<<24)){
			selectedKernel = dpKernels[6];
		}
		else if(p <=(1<<25)){
			selectedKernel = dpKernels[5];
		}
		else{
			selectedKernel = dpKernels[4];
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
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

		return D;
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<float> Field
	 * D = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::muladd<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& D, const float& beta, const BlasMatrix<Modular<float> >& C,
		const float& alpha, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[4];
		kernelsAvailable &= spKernelsAvailable[5];
		kernelsAvailable &= spKernelsAvailable[6];
		kernelsAvailable &= spKernelsAvailable[7];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<float>, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >()(
				_F,D,beta,C,alpha,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			return BlasMatrixDomainMulAdd<
				Modular<float>, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >()(
				_F,D,beta,C,alpha,A,B);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(C);

		float p = _F.characteristic();
		float tempAlpha = fmod(alpha, p);
		float tempBeta = fmod(beta, p);

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[7];
		}
		else if(p <=(1<<9)){
			selectedKernel = spKernels[6];
		}
		else if(p <=(1<<10)){
			selectedKernel = spKernels[5];
		}
		else{
			selectedKernel = spKernels[4];
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
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::muladdin<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		const double& beta, BlasMatrix<Modular<double> >& C, const double& alpha,
		const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[4];
		kernelsAvailable &= dpKernelsAvailable[5];
		kernelsAvailable &= dpKernelsAvailable[6];
		kernelsAvailable &= dpKernelsAvailable[7];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<double>, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >()(
				_F,beta,C,alpha,A,B);
		}

		return muladd<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(C,beta,C,alpha,A,B);
	}

	/*
	 * Specialization of general matrix-matrix multiplication and
	 * addition with scaling over a Modular<float> Field
	 * Places the results into the first genreral dense matrix
	 * C = beta.C + alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::muladdin<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		const float& beta, BlasMatrix<Modular<float> >& C, const float& alpha,
		const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[4];
		kernelsAvailable &= spKernelsAvailable[5];
		kernelsAvailable &= spKernelsAvailable[6];
		kernelsAvailable &= spKernelsAvailable[7];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<float>, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >()(
				_F,beta,C,alpha,A,B);
		}

		return muladd<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(C,beta,C,alpha,A,B);
	}

	/*
	 * Specialization of multiplication with scaling over
	 * a Modular<double> Field
	 * C = alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::mul<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& C, const double& alpha, const BlasMatrix<Modular<double> >& A,
		const BlasMatrix<Modular<double> >& B) const{

		return muladdin<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(0,C,alpha,A,B);
	}

	/*
	 * Specialization of multiplication with scaling over
	 * a Modular<float> Field
	 * C = alpha.A*B
	 */
	template <>
	template <>
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::mul<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& C, const float& alpha, const BlasMatrix<Modular<float> >& A,
		const BlasMatrix<Modular<float> >& B) const{

		return muladdin<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(0,C,alpha,A,B);
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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::axpy<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& D, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B,
		const BlasMatrix<Modular<double> >& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[8];
		kernelsAvailable &= dpKernelsAvailable[9];
		kernelsAvailable &= dpKernelsAvailable[10];
		kernelsAvailable &= dpKernelsAvailable[11];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(
				_F,D,_One,C,_One,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			D = mul<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,A,B);
			return addin<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(C);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[11];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[10];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[9];
		}
		else{
			selectedKernel = dpKernels[8];
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
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::axpy<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& D, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B,
		const BlasMatrix<Modular<float> >& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[8];
		kernelsAvailable &= spKernelsAvailable[9];
		kernelsAvailable &= spKernelsAvailable[10];
		kernelsAvailable &= spKernelsAvailable[11];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(
				_F,D,_One,C,_One,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			D = mul<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,A,B);
			return addin<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(C);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[11];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[10];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[9];
		}
		else{
			selectedKernel = spKernels[8];
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
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::axpyin<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& C, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[8];
		kernelsAvailable &= dpKernelsAvailable[9];
		kernelsAvailable &= dpKernelsAvailable[10];
		kernelsAvailable &= dpKernelsAvailable[11];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return muladdin<BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >(_One,C,_One,A,B);
		}

		return axpy<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(C,A,B,C);
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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::axpyin<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& C, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[8];
		kernelsAvailable &= spKernelsAvailable[9];
		kernelsAvailable &= spKernelsAvailable[10];
		kernelsAvailable &= spKernelsAvailable[11];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return muladdin<BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >(_One,C,_One,A,B);
		}

		return axpy<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(C,A,B,C);
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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::maxpy<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& D, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B,
		const BlasMatrix<Modular<double> >& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[12];
		kernelsAvailable &= dpKernelsAvailable[13];
		kernelsAvailable &= dpKernelsAvailable[14];
		kernelsAvailable &= dpKernelsAvailable[15];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(
				_F,D,_One,C,_MOne,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			BlasMatrix<Modular<double> > T(_F,D.rowdim(),D.coldim());
			T = mul<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(T,A,B);
			return sub<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,C,T);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(C);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[15];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[14];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[13];
		}
		else{
			selectedKernel = dpKernels[12];
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
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::maxpy<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& D, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B,
		const BlasMatrix<Modular<float> >& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[12];
		kernelsAvailable &= spKernelsAvailable[13];
		kernelsAvailable &= spKernelsAvailable[14];
		kernelsAvailable &= spKernelsAvailable[15];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(
				_F,D,_One,C,_MOne,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			BlasMatrix<Modular<float> > T(_F,D.rowdim(),D.coldim());
			T = mul<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(T,A,B);
			return sub<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,C,T);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(C);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[15];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[14];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[13];
		}
		else{
			selectedKernel = spKernels[12];
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
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::maxpyin<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& C, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[12];
		kernelsAvailable &= dpKernelsAvailable[13];
		kernelsAvailable &= dpKernelsAvailable[14];
		kernelsAvailable &= dpKernelsAvailable[15];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return muladdin<BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >(_One,C,_MOne,A,B);
		}

		return maxpy<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(C,A,B,C);
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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::maxpyin<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& C, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[12];
		kernelsAvailable &= spKernelsAvailable[13];
		kernelsAvailable &= spKernelsAvailable[14];
		kernelsAvailable &= spKernelsAvailable[15];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return muladdin<BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >(_One,C,_MOne,A,B);
		}

		return maxpy<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(C,A,B,C);
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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::axmy<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& D, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B,
		const BlasMatrix<Modular<double> >& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_double, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[16];
		kernelsAvailable &= dpKernelsAvailable[17];
		kernelsAvailable &= dpKernelsAvailable[18];
		kernelsAvailable &= dpKernelsAvailable[19];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<double>,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >()(
				_F,D,_MOne,C,_One,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			D = mul<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,A,B);
			return subin<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(C);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^53
		if(p <= (1<<21)){
			selectedKernel = dpKernels[19];
		}
		else if(p <= (1<<24)){
			selectedKernel = dpKernels[18];
		}
		else if(p <= (1<<25)){
			selectedKernel = dpKernels[17];
		}
		else{
			selectedKernel = dpKernels[16];
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
		D = oclReadMatrixBuffer<cl_double, BlasMatrix<Modular<double> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::axmy<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& D, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B,
		const BlasMatrix<Modular<float> >& C) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<
			cl_float, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,A,B,C);

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[16];
		kernelsAvailable &= spKernelsAvailable[17];
		kernelsAvailable &= spKernelsAvailable[18];
		kernelsAvailable &= spKernelsAvailable[19];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return BlasMatrixDomainMulAdd<
				Modular<float>,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >()(
				_F,D,_MOne,C,_One,A,B);
		}
		//If the buffers use too much memory partition into blocks -- TODO
		if(!memLevelsAllowed){
			D = mul<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,A,B);
			return subin<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(D,C);
		}

		//Check dimensions
		linbox_check( A.coldim() == B.rowdim());
		linbox_check( C.rowdim() == A.rowdim());
		linbox_check( C.coldim() == B.coldim());
		linbox_check( D.rowdim() == C.rowdim());
		linbox_check( D.coldim() == C.coldim());

		//Lock the device
		pthread_mutex_lock(deviceLock);

		//Allocate buffers
		cl_mem bufferD = oclCreateMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(D);
		cl_mem bufferA = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(A);
		cl_mem bufferB = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(B);
		cl_mem bufferC = oclCreateAndLoadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(C);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
		//p^2 * n < 2^23
		if(p <= (1<<7)){
			selectedKernel = spKernels[19];
		}
		else if(p <= (1<<9)){
			selectedKernel = spKernels[18];
		}
		else if(p <= (1<<10)){
			selectedKernel = spKernels[17];
		}
		else{
			selectedKernel = spKernels[16];
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
		D = oclReadMatrixBuffer<cl_float, BlasMatrix<Modular<float> > >(bufferD, D);

		//Delete OpenCL buffers
		tempErrcode = clReleaseMemObject(bufferD);
		tempErrcode = clReleaseMemObject(bufferA);
		tempErrcode = clReleaseMemObject(bufferB);
		tempErrcode = clReleaseMemObject(bufferC);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Unlock the device
		pthread_mutex_unlock(deviceLock);

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
	BlasMatrix<Modular<double> >& OpenCLMatrixDomain<Modular<double> >::axmyin<
		BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(
		BlasMatrix<Modular<double> >& C, const BlasMatrix<Modular<double> >& A, const BlasMatrix<Modular<double> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = dpKernelsAvailable[16];
		kernelsAvailable &= dpKernelsAvailable[17];
		kernelsAvailable &= dpKernelsAvailable[18];
		kernelsAvailable &= dpKernelsAvailable[19];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !doubleSupported || !kernelsAvailable){
			return muladdin<BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> >,BlasMatrix<Modular<double> > >(_MOne,C,_One,A,B);
		}

		return axmy<BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> >, BlasMatrix<Modular<double> > >(C,A,B,C);
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
	BlasMatrix<Modular<float> >& OpenCLMatrixDomain<Modular<float> >::axmyin<
		BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(
		BlasMatrix<Modular<float> >& C, const BlasMatrix<Modular<float> >& A, const BlasMatrix<Modular<float> >& B) const{

		//Check if kernels are available
		bool kernelsAvailable = spKernelsAvailable[16];
		kernelsAvailable &= spKernelsAvailable[17];
		kernelsAvailable &= spKernelsAvailable[18];
		kernelsAvailable &= spKernelsAvailable[19];

		//If it is not capable or not setup properly use default implementation
		if(!setupCorrect || !kernelsAvailable){
			return muladdin<BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> >,BlasMatrix<Modular<float> > >(_MOne,C,_One,A,B);
		}

		return axmy<BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> >, BlasMatrix<Modular<float> > >(C,A,B,C);
	}

} //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_INL

