/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain-setup.inl
 * Copyright (C) 2012 Matthew Wezowicz
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

#ifndef __LINBOX_opencl_matrix_domain_util_INL
#define __LINBOX_opencl_matrix_domain_util_INL

#include "linbox/algorithms/opencl-domain-factory.h"

namespace LinBox{

	/**
	 * @internal
	 * Initializes the OpenCL compute environment
	 */
	template <class Field>
	void OpenCLMatrixDomain<Field>::oclMatrixDomainInit(){

		OpenCLMatrixDomainFactory::oclMatrixDomainCreate(this);
	}

	/**
	 * @internal
	 * Releases OpenCL cumpute resources
	 */
	template <class Field>
	void OpenCLMatrixDomain<Field>::oclMatrixDomainRelease(unsigned int IDnum){

		OpenCLMatrixDomainFactory::oclMatrixDomainDestroy(IDnum);
	}

	/**
	 * @internal
	 * Checks to see if the memory levels required are possible
	 */
	template <class Field>
	template <typename T, class Operand1, class Operand2, class Operand3>
	bool OpenCLMatrixDomain<Field>::oclMemCheck(Operand1 &C, const Operand2 &A, const Operand3 &B) const{

		//Calculate dimensions after padding of matrices
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newCDimX = ((C.coldim() + 15) / 16) * 16;
		int newCDimY = ((C.rowdim() + 15) / 16) * 16;
		int newADimX = ((A.coldim() + 15) / 16) * 16;
		int newADimY = ((A.rowdim() + 15) / 16) * 16;
		int newBDimX = ((B.coldim() + 15) / 16) * 16;
		int newBDimY = ((B.rowdim() + 15) / 16) * 16;

		//Determine if each individual matrix will fit in a buffer
		bool temp = (maxBufferSize >= (newCDimX * newCDimY * sizeof(T)));
		temp &= (maxBufferSize >= (newADimX) * newADimY * sizeof(T));
		temp &= (maxBufferSize >= (newBDimX * newBDimY) * sizeof(T));

		//Determine if all three buffers will fit at the same time
		temp &= (memCapacity >= ((newCDimX * newCDimY) + (newADimX * newADimY) +
			(newBDimX * newBDimY)) * sizeof(T));

		return temp;
	}

	template <class Field>
	template <typename T, class Operand1, class Operand2, class Operand3>
	bool OpenCLMatrixDomain<Field>::oclMemCheck(Operand1& D, const Operand2& A, const Operand3& B,
		const Operand1& C) const{

		//Calculate dimensions after padding of matrices
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDDimX = ((D.coldim() + 15) / 16) * 16;
		int newDDimY = ((D.rowdim() + 15) / 16) * 16;
		int newADimX = ((A.coldim() + 15) / 16) * 16;
		int newADimY = ((A.rowdim() + 15) / 16) * 16;
		int newBDimX = ((B.coldim() + 15) / 16) * 16;
		int newBDimY = ((B.rowdim() + 15) / 16) * 16;
		int newCDimX = ((C.coldim() + 15) / 16) * 16;
		int newCDimY = ((C.rowdim() + 15) / 16) * 16;

		//Determine if each individual matrix will fit in a buffer
		bool temp = (maxBufferSize >= (newDDimX * newDDimY * sizeof(T)));
		temp &= (maxBufferSize >= (newADimX) * newADimY * sizeof(T));
		temp &= (maxBufferSize >= (newBDimX * newBDimY) * sizeof(T));
		temp &= (maxBufferSize >= (newCDimX * newCDimY) * sizeof(T));

		//Determine if all three buffers will fit at the same time
		temp &= (memCapacity >= ((newDDimX * newDDimY) + (newADimX * newADimY) +
			(newBDimX * newBDimY) + (newCDimX * newCDimY)) * sizeof(T));

		return temp;
	}

	/**
	 * @internal
	 * Functions to call the passed kernel on the passed buffers
	 */
	template <class Field>
	template <typename T, typename U>
	void OpenCLMatrixDomain<Field>::oclCallKernel(cl_mem bufferC, cl_mem bufferA, cl_mem bufferB,
		int widthA, int heightA ,int widthB, T p, cl_kernel selectedKernel) const{

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(U), (void*)&p);
		////updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		////updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually
	}

	template <class Field>
	template <typename T, typename U>
	void OpenCLMatrixDomain<Field>::oclCallKernel(cl_mem bufferD, cl_mem bufferA, cl_mem bufferB,
		cl_mem bufferC, int widthA, int heightA, int widthB, T p, cl_kernel selectedKernel) const{

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(U), (void*)&p);
		////updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

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
		////updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually
	}

	template <class Field>
	template <typename T, typename U>
	void OpenCLMatrixDomain<Field>::oclCallKernel(cl_mem bufferD, cl_mem bufferA, cl_mem bufferB,
		cl_mem bufferC, T alpha, T beta, int widthA, int heightA, int widthB, T p, cl_kernel selectedKernel) const{

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(selectedKernel, 0, sizeof(cl_mem), (void*)&bufferD);
		tempErrcode = clSetKernelArg(selectedKernel, 1, sizeof(U), (void*)&alpha);
		tempErrcode = clSetKernelArg(selectedKernel, 2, sizeof(cl_mem), (void*)&bufferA);
		tempErrcode = clSetKernelArg(selectedKernel, 3, sizeof(cl_mem), (void*)&bufferB);
		tempErrcode = clSetKernelArg(selectedKernel, 4, sizeof(U), (void*)&beta);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(cl_mem), (void*)&bufferC);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(cl_int), (void*)&widthA);
		tempErrcode = clSetKernelArg(selectedKernel, 7, sizeof(cl_int), (void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 8, sizeof(U), (void*)&p);
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
		////updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually
	}


	/**
	 * @internal
	 * Functions to partition the matrices into submatrix views
	 */
	template <class Field>
	template <class Operand1, class Operand2, class Operand3>
	std::vector<int> OpenCLMatrixDomain<Field>::oclPartition(Operand1& C, const Operand2& A, const Operand3& B,
		SubmatrixVector& VC, SubmatrixVector& VA, SubmatrixVector& VB) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<Element,Operand1,Operand2,Operand3>(C,A,B,C);

		std::vector<int> temp;

		if(memLevelsAllowed){
			// Create Submatrix views
			BlasSubmatrix<Field> SC(C);
			BlasSubmatrix<Field> SA(A);
			BlasSubmatrix<Field> SB(B);

			// Place Submatrices at the beginning of the vectors
			VC.push_back(SC);
			VA.push_back(SA);
			VB.push_back(SB);

			// Return the block dimensions
			temp.push_back(1); //CBlocksX
			temp.push_back(1); // CBlocksY
			temp.push_back(1); //ABlocksY
			temp.push_back(1); //ABlocksY
			temp.push_back(1); //CBlocksY
			temp.push_back(1); //CBlocksY

			return temp;
		}

		return temp;
	}

	template <class Field>
	template <class Operand1, class Operand2, class Operand3>
	std::vector<int> OpenCLMatrixDomain<Field>::oclPartition(Operand1& D, const Operand2& A, const Operand3& B,
		const Operand1& C, SubmatrixVector& VD, SubmatrixVector& VA, SubmatrixVector& VB, SubmatrixVector& VC) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<Element,Operand1,Operand2,Operand3>(D,A,B,C);

		std::vector<int> temp;

		if(memLevelsAllowed){
			// Create Submatrix views
			BlasSubmatrix<Field> SD(D);
			BlasSubmatrix<Field> SA(A);
			BlasSubmatrix<Field> SB(B);
			BlasSubmatrix<Field> SC(C);

			// Place Submatrices at the beginning of the vectors
			VD.push_back(SD);
			VA.push_back(SA);
			VB.push_back(SB);
			VC.push_back(SC);

			// Return the block dimensions
			temp.push_back(1); //DBlocksX & CBlocksX
			temp.push_back(1); //DBlocksY & CBlocksY
			temp.push_back(1); //ABlocksY
			temp.push_back(1); //ABlocksY
			temp.push_back(1); //CBlocksY
			temp.push_back(1); //CBlocksY

			return temp;
		}

		return temp;
	}

	//Prints the appropriate error message
	template <class Field>
	void OpenCLMatrixDomain<Field>::printClErrstring(cl_int err) const{
		switch (err) {
		case CL_SUCCESS:
			std::cout << "Success!\n";
			break;
		case CL_DEVICE_NOT_FOUND:
			std::cout << "Device not found.\n";
			break;
		case CL_DEVICE_NOT_AVAILABLE:
			std::cout << "Device not available\n";
			break;
		case CL_COMPILER_NOT_AVAILABLE:
			std::cout << "Compiler not available\n";
			break;
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:
			std::cout << "Memory object allocation failure\n";
			break;
		case CL_OUT_OF_RESOURCES:
			std::cout << "Out of resources\n";
			break;
		case CL_OUT_OF_HOST_MEMORY:
			std::cout << "Out of host memory\n";
			break;
		case CL_PROFILING_INFO_NOT_AVAILABLE:
			std::cout << "Profiling information not available\n";
			break;
		case CL_MEM_COPY_OVERLAP:
			std::cout << "Memory copy overlap\n";
			break;
		case CL_IMAGE_FORMAT_MISMATCH:
			std::cout << "Image format mismatch\n";
			break;
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:
			std::cout << "Image format not supported\n";
			break;
		case CL_BUILD_PROGRAM_FAILURE:
			std::cout << "Program build failure\n";
			break;
		case CL_MAP_FAILURE:
			std::cout << "Map failure\n";
			break;
		case CL_INVALID_VALUE:
			std::cout << "Invalid value\n";
			break;
		case CL_INVALID_DEVICE_TYPE:
			std::cout << "Invalid device type\n";
			break;
		case CL_INVALID_PLATFORM:
			std::cout << "Invalid platform\n";
			break;
		case CL_INVALID_DEVICE:
			std::cout << "Invalid device\n";
			break;
		case CL_INVALID_CONTEXT:
			std::cout << "Invalid context\n";
			break;
		case CL_INVALID_QUEUE_PROPERTIES:
			std::cout << "Invalid queue properties\n";
			break;
		case CL_INVALID_COMMAND_QUEUE:
			std::cout << "Invalid command queue\n";
			break;
		case CL_INVALID_HOST_PTR:
			std::cout << "Invalid host pointer\n";
			break;
		case CL_INVALID_MEM_OBJECT:
			std::cout << "Invalid memory object\n";
			break;
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:
			std::cout << "Invalid image format descriptor\n";
			break;
		case CL_INVALID_IMAGE_SIZE:
			std::cout << "Invalid image size\n";
			break;
		case CL_INVALID_SAMPLER:
			std::cout << "Invalid sampler\n";
			break;
		case CL_INVALID_BINARY:
			std::cout << "Invalid binary\n";
			break;
		case CL_INVALID_BUILD_OPTIONS:
			std::cout << "Invalid build options\n";
			break;
		case CL_INVALID_PROGRAM:
			std::cout << "Invalid program\n";
			break;
		case CL_INVALID_PROGRAM_EXECUTABLE:
			std::cout << "Invalid program executable\n";
			break;
		case CL_INVALID_KERNEL_NAME:
			std::cout << "Invalid kernel name\n";
			break;
		case CL_INVALID_KERNEL_DEFINITION:
			std::cout << "Invalid kernel definition\n";
			break;
		case CL_INVALID_KERNEL:
			std::cout << "Invalid kernel\n";
			break;
		case CL_INVALID_ARG_INDEX:
			std::cout << "Invalid argument index\n";
			break;
		case CL_INVALID_ARG_VALUE:
			std::cout << "Invalid argument value\n";
			break;
		case CL_INVALID_ARG_SIZE:
			std::cout << "Invalid argument size\n";
			break;
		case CL_INVALID_KERNEL_ARGS:
			std::cout << "Invalid kernel arguments\n";
			break;
		case CL_INVALID_WORK_DIMENSION:
			std::cout << "Invalid work dimension\n";
			break;
		case CL_INVALID_WORK_GROUP_SIZE:
			std::cout << "Invalid work group size\n";
			break;
		case CL_INVALID_WORK_ITEM_SIZE:
			std::cout << "Invalid work item size\n";
			break;
		case CL_INVALID_GLOBAL_OFFSET:
			std::cout << "Invalid global offset\n";
			break;
		case CL_INVALID_EVENT_WAIT_LIST:
			std::cout << "Invalid event wait list\n";
			break;
		case CL_INVALID_EVENT:
			std::cout << "Invalid event\n";
			break;
		case CL_INVALID_OPERATION:
			std::cout << "Invalid operation\n";
			break;
		case CL_INVALID_GL_OBJECT:
			std::cout << "Invalid OpenGL object\n";
			break;
		case CL_INVALID_BUFFER_SIZE:
			std::cout << "Invalid buffer size\n";
			break;
		case CL_INVALID_MIP_LEVEL:
			std::cout << "Invalid mip-map level\n";
			break;
		default:
			std::cout << "Unknown\n";
			break;
		}
	}

} //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_util_INL