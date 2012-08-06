/* linbox/algorithms/opencl-domain-setup.inl
 * Copyright (C) 2012 Matthew Wezowicz
 *
 * Written by Matthew Wezowicz <mwezz@udel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_opencl_matrix_domain_util_INL
#define __LINBOX_opencl_matrix_domain_util_INL

#include <cstdio>
#include <utility>
#include "linbox/algorithms/opencl-domain-factory.h"

namespace LinBox{

	/**
	 * @internal
	 * Initializes the OpenCL compute environment
	 */
	template <class Field>
	void OpenCLMatrixDomain<Field>::oclMatrixDomainAcquire(){

		OpenCLMatrixDomainFactory::oclMatrixDomainInstance(this);
	}

	/**
	 * @internal
	 * Releases OpenCL cumpute resources
	 */
	template <class Field>
	void OpenCLMatrixDomain<Field>::oclMatrixDomainRelease(unsigned int IDnum){

		OpenCLMatrixDomainFactory::oclMatrixDomainDeallocate(IDnum);
	}

	/**
	 * @internal
	 * Checks to see if the memory levels required are possible
	 */
	template <>
	template <class Operand1, class Operand2, class Operand3>
	bool OpenCLMatrixDomain<Modular<double> >::oclMemCheck(
		Operand1& D,
		const Operand2& A,
		const Operand3& B,
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
		bool temp = (maxBufferSize >= (newDDimX * newDDimY * sizeof(cl_double)));
		temp &= (maxBufferSize >= (newADimX) * newADimY * sizeof(cl_double));
		temp &= (maxBufferSize >= (newBDimX * newBDimY) * sizeof(cl_double));
		temp &= (maxBufferSize >= (newCDimX * newCDimY) * sizeof(cl_double));

		//Determine if all three buffers will fit at the same time
		temp &= (memCapacity >= ((newDDimX * newDDimY) + (newADimX * newADimY) +
			      (newBDimX * newBDimY) + (newCDimX * newCDimY)) * sizeof(cl_double));

		return temp;
	}

	template <>
	template <class Operand1, class Operand2, class Operand3>
	bool OpenCLMatrixDomain<Modular<float> >::oclMemCheck(
		Operand1& D,
		const Operand2& A,
		const Operand3& B,
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
		bool temp = (maxBufferSize >= (newDDimX * newDDimY * sizeof(cl_float)));
		temp &= (maxBufferSize >= (newADimX) * newADimY * sizeof(cl_float));
		temp &= (maxBufferSize >= (newBDimX * newBDimY) * sizeof(cl_float));
		temp &= (maxBufferSize >= (newCDimX * newCDimY) * sizeof(cl_float));

		//Determine if all three buffers will fit at the same time
		temp &= (memCapacity >= ((newDDimX * newDDimY) + (newADimX * newADimY) +
			      (newBDimX * newBDimY) + (newCDimX * newCDimY)) * sizeof(cl_float));

		return temp;
	}

	template <>
	template <>
	bool OpenCLMatrixDomain<Modular<double> >::oclMemCheck<std::pair<int,int> >(
		std::pair<int,int>& D,
		std::pair<int,int>& A,
		std::pair<int,int>& B,
		std::pair<int,int>& C) const{

		//Calculate dimensions after padding of matrices
		//((A.second / 16) + (A.second % 16 == 0 ? 0 : 1)) * 16
		int newDDimX = ((D.second + 15) / 16) * 16;
		int newDDimY = ((D.first + 15) / 16) * 16;
		int newADimX = ((A.second + 15) / 16) * 16;
		int newADimY = ((A.first + 15) / 16) * 16;
		int newBDimX = ((B.second + 15) / 16) * 16;
		int newBDimY = ((B.first + 15) / 16) * 16;
		int newCDimX = ((C.second + 15) / 16) * 16;
		int newCDimY = ((C.first + 15) / 16) * 16;

		//Determine if each individual matrix will fit in a buffer
		bool temp = (maxBufferSize >= (newDDimX * newDDimY * sizeof(cl_double)));
		temp &= (maxBufferSize >= (newADimX) * newADimY * sizeof(cl_double));
		temp &= (maxBufferSize >= (newBDimX * newBDimY) * sizeof(cl_double));
		temp &= (maxBufferSize >= (newCDimX * newCDimY) * sizeof(cl_double));

		//Determine if all three buffers will fit at the same time
		temp &= (memCapacity >= ((newDDimX * newDDimY) + (newADimX * newADimY) +
			      (newBDimX * newBDimY) + (newCDimX * newCDimY)) * sizeof(cl_double));

		return temp;
	}

	template <>
	template <>
	bool OpenCLMatrixDomain<Modular<float> >::oclMemCheck<std::pair<int,int> >(
		std::pair<int,int>& D,
		std::pair<int,int>& A,
		std::pair<int,int>& B,
		std::pair<int,int>& C) const{

		//Calculate dimensions after padding of matrices
		//((A.second / 16) + (A.second % 16 == 0 ? 0 : 1)) * 16
		int newDDimX = ((D.second + 15) / 16) * 16;
		int newDDimY = ((D.first + 15) / 16) * 16;
		int newADimX = ((A.second + 15) / 16) * 16;
		int newADimY = ((A.first + 15) / 16) * 16;
		int newBDimX = ((B.second + 15) / 16) * 16;
		int newBDimY = ((B.first + 15) / 16) * 16;
		int newCDimX = ((C.second + 15) / 16) * 16;
		int newCDimY = ((C.first + 15) / 16) * 16;

		//Determine if each individual matrix will fit in a buffer
		bool temp = (maxBufferSize >= (newDDimX * newDDimY * sizeof(cl_float)));
		temp &= (maxBufferSize >= (newADimX) * newADimY * sizeof(cl_float));
		temp &= (maxBufferSize >= (newBDimX * newBDimY) * sizeof(cl_float));
		temp &= (maxBufferSize >= (newCDimX * newCDimY) * sizeof(cl_float));

		//Determine if all three buffers will fit at the same time
		temp &= (memCapacity >= ((newDDimX * newDDimY) + (newADimX * newADimY) +
			      (newBDimX * newBDimY) + (newCDimX * newCDimY)) * sizeof(cl_float));

		return temp;
	}

	/**
	 * @internal
	 * Functions to call the passed kernel on the passed buffers
	 */
	template <class Field>
	template <typename T, typename U>
	void OpenCLMatrixDomain<Field>::oclCallKernel(
		cl_mem bufferC,
		cl_mem bufferA,
		cl_mem bufferB,
		int widthA,
		int heightA,
		int widthB,
		T p,
		cl_kernel selectedKernel) const{

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(
			selectedKernel,
			0,
			sizeof(cl_mem),
			(void*)&bufferC);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			1,
			sizeof(cl_mem),
			(void*)&bufferA);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			2,
			sizeof(cl_mem),
			(void*)&bufferB);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			3,
			sizeof(cl_int),
			(void*)&widthA);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			4,
			sizeof(cl_int),
			(void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 5, sizeof(U), (void*)&p);

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(
			commandQue,
			selectedKernel,
			2,
			NULL,
			globalWorkSize,
			localWorkSize,
			0,
			NULL,
			NULL);
	}

	template <class Field>
	template <typename T, typename U>
	void OpenCLMatrixDomain<Field>::oclCallKernel(
		cl_mem bufferD,
		cl_mem bufferA,
		cl_mem bufferB,
		cl_mem bufferC,
		int widthA,
		int heightA,
		int widthB,
		T p,
		cl_kernel selectedKernel) const{

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(
			selectedKernel,
			0,
			sizeof(cl_mem),
			(void*)&bufferD);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			1,
			sizeof(cl_mem),
			(void*)&bufferA);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			2,
			sizeof(cl_mem),
			(void*)&bufferB);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			3,
			sizeof(cl_mem),
			(void*)&bufferC);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			4,
			sizeof(cl_int),
			(void*)&widthA);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			5,
			sizeof(cl_int),
			(void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 6, sizeof(U), (void*)&p);

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(
			commandQue,
			selectedKernel,
			2,
			NULL,
			globalWorkSize,
			localWorkSize,
			0,
			NULL,
			NULL);
	}

	template <class Field>
	template <typename T, typename U>
	void OpenCLMatrixDomain<Field>::oclCallKernel(
		cl_mem bufferD,
		cl_mem bufferA,
		cl_mem bufferB,
		cl_mem bufferC,
		T alpha,
		T beta,
		int widthA,
		int heightA,
		int widthB,
		T p,
		cl_kernel selectedKernel) const{

		//Pass kernel arguments
		cl_int tempErrcode;
		tempErrcode = clSetKernelArg(
			selectedKernel,
			0,
			sizeof(cl_mem),
			(void*)&bufferD);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			1,
			sizeof(U),
			(void*)&alpha);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			2,
			sizeof(cl_mem),
			(void*)&bufferA);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			3,
			sizeof(cl_mem),
			(void*)&bufferB);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			4,
			sizeof(U),
			(void*)&beta);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			5,
			sizeof(cl_mem),
			(void*)&bufferC);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			6,
			sizeof(cl_int),
			(void*)&widthA);
		tempErrcode = clSetKernelArg(
			selectedKernel,
			7,
			sizeof(cl_int),
			(void*)&widthB);
		tempErrcode = clSetKernelArg(selectedKernel, 8, sizeof(U), (void*)&p);

		//Set NDRange
		size_t localWorkSize[2];
		size_t globalWorkSize[2];
		localWorkSize[0] = 16;
		localWorkSize[1] = 16;
		globalWorkSize[0] = widthB;
		globalWorkSize[1] = heightA;

		//Launch kernel
		tempErrcode = clEnqueueNDRangeKernel(
			commandQue,
			selectedKernel,
			2,
			NULL,
			globalWorkSize,
			localWorkSize,
			0,
			NULL,
			NULL);
	}


	/**
	 * @internal
	 * Functions to partition the matrices into submatrix views
	 */
	template <class Field>
	template <class Operand1, class Operand2, class Operand3>
	std::vector<int> OpenCLMatrixDomain<Field>::oclPartition(
		Operand1& C,
		const Operand2& A,
		const Operand3& B,
		std::vector<SubmatrixAdapter<Operand1> >& VC,
		std::vector<SubmatrixAdapter<Operand2> >& VA,
		std::vector<SubmatrixAdapter<Operand3> >& VB) const{

		std::vector<SubmatrixAdapter<Operand1> > VT;
		return oclPartition<Operand1,Operand2,Operand3>(C,A,B,C,VT,VA,VB,VC);
	}

	template <class Field>
	template <class Operand1, class Operand2, class Operand3>
	std::vector<int> OpenCLMatrixDomain<Field>::oclPartition(
		Operand1& D,
		const Operand2& A,
		const Operand3& B,
		const Operand1& C,
		std::vector<SubmatrixAdapter<Operand1> >& VD,
		std::vector<SubmatrixAdapter<Operand2> >& VA,
		std::vector<SubmatrixAdapter<Operand3> >& VB,
		std::vector<SubmatrixAdapter<Operand1> >& VC) const{

		//Compute if the OpenCL device is capable of working with the required ammounts of memory
		bool memLevelsAllowed = oclMemCheck<Operand1,Operand2,Operand3>(D,A,B,C);

		std::vector<int> temp;

		if(memLevelsAllowed){
			// Create Submatrix views
			SubmatrixAdapter<Operand1> SD(D);
			SubmatrixAdapter<Operand2> SA(A);
			SubmatrixAdapter<Operand3> SB(B);
			SubmatrixAdapter<Operand1> SC(C);

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

		//Begin search for Submatrices small enough to search on the device
		int divisionFactor = 1;
		int DRows = D.rowdim();
		int DCols = D.coldim();
		int ARows = A.rowdim();
		int ACols = A.rowdim();
		int BRows = B.rowdim();
		int BCols = B.coldim();

		if(true){ //Default partitioning scheme

			//Loop until Submatrices that fit on the device have been found
			while(!memLevelsAllowed){
				//Increase the number of subsections
				divisionFactor++;

				//Compute the number of rows and columns in each subsections
				int subRows = DRows / divisionFactor;
				int subCols = DCols / divisionFactor;
				int subSharedDim = ACols / divisionFactor;

				//Check if adjustment is need for some of the subsections
				bool addToSubRows = false;
				bool addToSubCols = false;
				bool addToSubSharedDim = false;

				if((subRows * divisionFactor) != DRows){
					addToSubRows = true;
				}
				if((subCols * divisionFactor) != DCols){
					addToSubCols = true;
				}
				if((subSharedDim * divisionFactor) != ACols){
					addToSubSharedDim = true;
				}

				//Determine of the largest subsections will fit on the device together
				int largestSubRows = (addToSubRows ? subRows : (subRows + 1));
				int largestSubCols = (addToSubCols ? subCols : (subCols + 1));
				int largestSubSharedDim = (addToSubSharedDim ? subSharedDim :
				                                               (subSharedDim + 1));

				std::pair<int,int> largestSubD(largestSubRows, largestSubCols);
				std::pair<int,int> largestSubA(largestSubRows, largestSubSharedDim);
				std::pair<int,int> largestSubB(largestSubSharedDim, largestSubCols);
				std::pair<int,int> largestSubC(largestSubRows, largestSubCols);

				memLevelsAllowed = oclMemCheck<std::pair<int,int> >(
					largestSubD,
					largestSubA,
					largestSubB,
					largestSubC);

				//If the largest subsections can fit on the device together
				//Begin partitioning the input matrices
				if(memLevelsAllowed){
					// Return the block dimensions
					temp.push_back(divisionFactor); //DBlocksX & CBlocksX
					temp.push_back(divisionFactor); //DBlocksY & CBlocksY
					temp.push_back(divisionFactor); //ABlocksY
					temp.push_back(divisionFactor); //ABlocksY
					temp.push_back(divisionFactor); //CBlocksY
					temp.push_back(divisionFactor); //CBlocksY

					//Loop through all but the last row
					for(int blockY = 0; blockY < (divisionFactor - 1); blockY++){
						int DRowsOffset = blockY * subRows;
						int ARowsOffset = blockY * subRows;
						int BRowsOffset = blockY * subSharedDim;

						//Loop through all but the last column
						for(int blockX = 0; blockX < (divisionFactor - 1); blockX++){
							int DColsOffset = blockX * subCols;
							int AColsOffset = blockX * subSharedDim;
							int BColsOffset = blockX * subCols;

							SubmatrixAdapter<Operand1> SD(
								D,
								DRowsOffset,
								DColsOffset,
								subRows,
								subCols);
							SubmatrixAdapter<Operand2> SA(
								A,
								ARowsOffset,
								AColsOffset,
								subRows,
								subSharedDim);
							SubmatrixAdapter<Operand3> SB(
								B,
								BRowsOffset,
								BColsOffset,
								subSharedDim,
								subCols);
							SubmatrixAdapter<Operand1> SC(
								C,
								DRowsOffset,
								DColsOffset,
								subRows,
								subCols);

							VD.push_back(SD);
							VA.push_back(SA);
							VB.push_back(SB);
							VC.push_back(SC);
						}

						int DColsOffset = (divisionFactor - 1) * subCols;
						int AColsOffset = (divisionFactor - 1) * subSharedDim;
						int BColsOffset = (divisionFactor - 1) * subCols;

						SubmatrixAdapter<Operand1> SD(
							D,
							DRowsOffset,
							DColsOffset,
							subRows,
							(DCols - DColsOffset));
						SubmatrixAdapter<Operand2> SA(
							A,
							ARowsOffset,
							AColsOffset,
							subRows,
							(ACols - AColsOffset));
						SubmatrixAdapter<Operand3> SB(
							B,
							BRowsOffset,
							BColsOffset,
							subSharedDim,
							(BCols - BColsOffset));
						SubmatrixAdapter<Operand1> SC(
							C,
							DRowsOffset,
							DColsOffset,
							subRows,
							(DCols - DColsOffset));

						VD.push_back(SD);
						VA.push_back(SA);
						VB.push_back(SB);
						VC.push_back(SC);
					}

					//Partition the last row
					int DRowsOffset = (divisionFactor - 1) * subRows;
					int ARowsOffset = (divisionFactor - 1) * subRows;
					int BRowsOffset = (divisionFactor - 1) * subSharedDim;

					for(int blockX = 0; blockX < (divisionFactor - 1); blockX++){
						int DColsOffset = blockX * subCols;
						int AColsOffset = blockX * subSharedDim;
						int BColsOffset = blockX * subCols;

						SubmatrixAdapter<Operand1> SD(
							D,
							DRowsOffset,
							DColsOffset,
							(DRows - DRowsOffset),
							subCols);
						SubmatrixAdapter<Operand2> SA(
							A,
							ARowsOffset,
							AColsOffset,
							(ARows - ARowsOffset),
							subSharedDim);
						SubmatrixAdapter<Operand3> SB(
							B,
							BRowsOffset,
							BColsOffset,
							(BRows - BRowsOffset),
							subCols);
						SubmatrixAdapter<Operand1> SC(
							C,
							DRowsOffset,
							DColsOffset,
							(DRows - DRowsOffset),
							subCols);

						VD.push_back(SD);
						VA.push_back(SA);
						VB.push_back(SB);
						VC.push_back(SC);
					}

					int DColsOffset = (divisionFactor - 1) * subCols;
					int AColsOffset = (divisionFactor - 1) * subSharedDim;
					int BColsOffset = (divisionFactor - 1) * subCols;

					SubmatrixAdapter<Operand1> SD(
						D,
						DRowsOffset,
						DColsOffset,
						(DRows - DRowsOffset),
						(DCols - DColsOffset));
					SubmatrixAdapter<Operand2> SA(
						A,
						ARowsOffset,
						AColsOffset,
						(ARows - ARowsOffset),
						(ACols - AColsOffset));
					SubmatrixAdapter<Operand3> SB(
						B,
						BRowsOffset,
						BColsOffset,
						(BRows - BRowsOffset),
						(BCols - BColsOffset));
					SubmatrixAdapter<Operand1> SC(
						C,
						DRowsOffset,
						DColsOffset,
						(DRows - DRowsOffset),
						(DCols - DColsOffset));

					VD.push_back(SD);
					VA.push_back(SA);
					VB.push_back(SB);
					VC.push_back(SC);

				}
			}
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
			std::cout << "Device not found\n";
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

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
