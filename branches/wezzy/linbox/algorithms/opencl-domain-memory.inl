/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/opencl-domain-setup.inl
 * Copyright (C) 2011 Matthew Wezowicz
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

#ifndef __LINBOX_opencl_matrix_domain_memory_INL
#define __LINBOX_opencl_matrix_domain_memory_INL

#include <new>

#include "CL/cl.hpp"

namespace LinBox{

	/**
	 * @internal
	 * Checks to see if the memory levels required are possible
	 */
	template<class Field>
	template<typename T, class Operand1, class Operand2, class Operand3>
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

	template<class Field>
	template<typename T, class Operand1, class Operand2, class Operand3>
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
	 * Pads a BlasMatrix into a form appropriate for OpenCL use and returns the OpenCL buffer
	 */
	template<class Field>
	template<typename T, class Operand1>
	cl_mem OpenCLMatrixDomain<Field>::oclPadMatrix(cl_mem matrixBuffer, int matrixBufferSize,
		int newDimX, const Operand1 &matrix) const{

		//Set starting positions
		int matrixBufferPosition = 0;
		int dataOffset = 0;

		//Allocates a 32mb buffer for padding
		T* paddingBuffer = (T*)operator new(32 * 1024 * 1024);

		//Calculate the size of the padding buffer in number of elements
		const int paddingBufferSize = (32 * 1024 * 1024 / sizeof(T));

		//Loops while there is still space in the matrixBuffer
		while(matrixBufferPosition < matrixBufferSize){

			//Sets the starting position in the buffer
			int paddingBufferPosition = 0;

			//Loops while there is still space in the buffer
			while(paddingBufferPosition < paddingBufferSize){
				int count = 0;

				//Puts one row of data into the buffer while there is space in the buffer and data left
				while(count < (int)matrix.coldim() && paddingBufferPosition < paddingBufferSize
					&& dataOffset < (int)(matrix.coldim() * matrix.rowdim())){

					//Put entry of matrix into the padding buffer
					paddingBuffer[paddingBufferPosition] = matrix.getEntry(
						(dataOffset / matrix.coldim()),(dataOffset % matrix.coldim()));

					//Increment the count for the row, paddingBuffer, and matrix
					count++;
					paddingBufferPosition++;
					dataOffset++;
				}

				//Padds 0's until end of padded row while there is still space in the buffer
				while(count < newDimX && paddingBufferPosition < paddingBufferSize){

					//Place a zero into the paddingBuffer
					paddingBuffer[paddingBufferPosition] = 0.0;

					//Increment the count for the rwo and the paddingBuffer
					count++;
					paddingBufferPosition++;
				}
			}

			int transferSize;

			//Transfer the paddingBuffer to the matrixBuffer
			if((matrixBufferPosition + paddingBufferSize) <= matrixBufferSize){
				transferSize = (32 * 1024 * 1024);
			}
			//Transfer the partial paddingBuffer to the matrixBuffer
			else{
				transferSize = (matrixBufferSize - matrixBufferPosition) * sizeof(T);
			}

			cl_int tempErrcode;
			tempErrcode = clEnqueueWriteBuffer(commandQue, matrixBuffer, CL_TRUE,
				(matrixBufferPosition * sizeof(T)), transferSize, paddingBuffer, 0, NULL, NULL);
			//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

			//Increment position in matrixBuffer by the size of the paddingBuffer
			matrixBufferPosition += paddingBufferSize;
		}

		//Deallocates buffer
		delete paddingBuffer;

		return matrixBuffer;
	}

	/**
	 * @internal
	 * Dedads an OpenCL buffer and puts it back into a BlasMatrix
	 */
	template<class Field>
	template<typename T, class Operand1>
	Operand1& OpenCLMatrixDomain<Field>::oclDepadMatrix(cl_mem matrixBuffer, int matrixBufferSize,
		int outputSize, int newDimX, Operand1& matrix) const{

		//Set starting positions
		int matrixBufferPosition = 0;
		int dataOffset = 0;

		//Allocates a 32mb buffer for depadding
		T* depaddingBuffer = (T*)operator new(32 * 1024 * 1024);

		//Calculate the size of the depadding buffer in number of elements
		const int depaddingBufferSize = (32 * 1024 * 1024 / sizeof(T));

		//Loops while there are still elements in the matrixBuffer
		while(dataOffset < outputSize){
			int transferSize;

			//Transfer a full depaddingBuffer worth of elements back to the host
			if((matrixBufferPosition + depaddingBufferSize) <= matrixBufferSize){
				transferSize = (32 * 1024 * 1024);
			}
			//Transfer a partial depaddingBuffer worth of elements back to the host
			else{
				transferSize = (matrixBufferSize - matrixBufferPosition) * sizeof(T);
			}

			cl_int tempErrcode;
			tempErrcode = clEnqueueReadBuffer(commandQue, matrixBuffer, CL_TRUE,
				(matrixBufferPosition * sizeof(T)), transferSize, depaddingBuffer,
				0, NULL, NULL);
			//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

			//Set depaddiingBuffer start position
			int depaddingBufferPosition = 0;

			//Loops while there are still elements in the depaddingBuffer
			while(depaddingBufferPosition < depaddingBufferSize){
				int count = 0;

				//Puts one row of data into the matrix while there are elements in the depaddingBuffer
				while(count < (int)matrix.coldim() && depaddingBufferPosition < depaddingBufferSize
					&& dataOffset < outputSize){

					//Put entry of depadding buffer into the matrix
					matrix.setEntry((dataOffset / matrix.coldim()),(dataOffset % matrix.coldim()),
						depaddingBuffer[depaddingBufferPosition]);

					//Increment the count for the row, depadding buffer, and matrix
					count++;
					depaddingBufferPosition++;
					dataOffset++;
				}

				//Skip over the padding zero's
				while(count < newDimX && depaddingBufferPosition < depaddingBufferSize){
					count++;
					depaddingBufferPosition++;
				}
			}

			//Increment position in matrixBuffer by the size of the paddingBuffer
			matrixBufferPosition += depaddingBufferSize;
		}

		//Deallocates buffer
		delete depaddingBuffer;

		return matrix;
	}

	/**
	 * @internal
	 * Creates an empty matrix buffer on the OpenCL device from the dimensions of the matrix
	 */
	template<class Field>
	template<typename T, class Opernad1>
	cl_mem OpenCLMatrixDomain<Field>::oclCreateMatrixBuffer(Opernad1& matrix) const{

		//Calculate dimensions after padding of matrix
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDimX = ((matrix.coldim() + 15) / 16) * 16;
		int newDimY = ((matrix.rowdim() + 15) / 16) * 16;

		//Allocate buffer
		cl_int tempErrcode;
		cl_mem matrixBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE,
			(newDimX * newDimY * sizeof(T)), 0, &tempErrcode);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		return matrixBuffer;
	}

	/**
	 * @internal
	 * Creates a matrix buffer on the OpenCL device from the dimensions of the passed in matrix
	 * and load the contents of the matrix into the matrix buffer
	 */
	template<class Field>
	template<typename T, class Operand1>
	cl_mem OpenCLMatrixDomain<Field>::oclCreateAndLoadMatrixBuffer(const Operand1 &matrix) const{

		//Calculate dimensions after padding of matrix
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDimX = ((matrix.coldim() + 15) / 16) * 16;
		int newDimY = ((matrix.rowdim() + 15) / 16) * 16;

		//Allocate buffer
		cl_int tempErrcode;
		cl_mem matrixBuffer = clCreateBuffer(context, CL_MEM_READ_WRITE,
			(newDimX * newDimY * sizeof(T)), 0, &tempErrcode);
		//updateErrcode(tempErrcode); //Does not work because of const -- will fix eventually

		//Calculate number of elements in the matrixBuffer
		int matrixBufferSize = newDimX * newDimY;

		return oclPadMatrix<T, Operand1>(matrixBuffer, matrixBufferSize, newDimX, matrix);
	}

	/**
	 * @internal
	 * Read back the contents of the matrix buffer into the matrix and return a refence to the matrix
	 */
	template<class Field>
	template<typename T, class Operand2>
	Operand2& OpenCLMatrixDomain<Field>::oclReadMatrixBuffer(cl_mem matrixBuffer, Operand2 &matrix) const{

		//Calculate dimensions after padding of matrix
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDimX = ((matrix.coldim() + 15) / 16) * 16;
		int newDimY = ((matrix.rowdim() + 15) / 16) * 16;

		//Calculate number of elements in the matrixBuffer
		int matrixBufferSize = newDimX * newDimY;

		//Calculate number of elements in the matrix
		int outputSize = matrix.coldim() * matrix.rowdim();

		return oclDepadMatrix<T, Operand2>(matrixBuffer, matrixBufferSize, outputSize,
			newDimX, matrix);
	}

}; //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_memory_INL