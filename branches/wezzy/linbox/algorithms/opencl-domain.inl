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
	 * multiplying two general dense matrices.
	 * over a Modular<double> Field.
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

		cl_mem bufferC = createMatrixBuffer<BlasMatrix<double> >(C);
		cl_mem bufferA = createAndLoadMatrixBuffer<BlasMatrix<double> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<BlasMatrix<double> >(B);

		double p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
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

		C = readMatrixBuffer<BlasMatrix<double> >(bufferC, C);

		return C;
	}

	/* Specialization of Mul for
	 * multiplying two general dense matrices.
	 * over a Modular<float> Field.
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

		cl_mem bufferC = createMatrixBuffer<BlasMatrix<float> >(C);
		cl_mem bufferA = createAndLoadMatrixBuffer<BlasMatrix<float> >(A);
		cl_mem bufferB = createAndLoadMatrixBuffer<BlasMatrix<float> >(B);

		float p = _F.characteristic();

		cl_kernel selectedKernel;

		// Select OpenCL kernel based on the size of the modulus factor for maximum performance
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

		C = readMatrixBuffer<BlasMatrix<float> >(bufferC, C);

		return C;
	}

} //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_INL

