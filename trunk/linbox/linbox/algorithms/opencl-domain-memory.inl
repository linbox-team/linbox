/* linbox/algorithms/opencl-domain-setup.inl
 * Copyright (C) 2011-2012 Matthew Wezowicz
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

#ifndef __LINBOX_opencl_matrix_domain_memory_INL
#define __LINBOX_opencl_matrix_domain_memory_INL

#include <cstdlib>

#include "CL/cl.h"

namespace LinBox{

	/**
	 * @internal
	 * Pads a BlasMatrix into a form appropriate for OpenCL use
	 * and returns the OpenCL buffer
	 */
	template <class Field>
	template <typename T, class Operand1>
	cl_mem OpenCLMatrixDomain<Field>::oclPadMatrix(
		cl_mem matrixBuffer,
		int matrixBufferSize,
		int newDimX,
		const Operand1 &matrix) const{

		//Set starting positions
		int matrixBufferPosition = 0;
		int dataOffset = 0;
		int elementsPadded = 0;

		//Allocates a 32mb buffer for padding
		T* paddingBuffer = (T*)malloc(32 * 1024 * 1024);

		//Calculate the size of the padding buffer in number of elements
		const int paddingBufferSize = (32 * 1024 * 1024 / (int)sizeof(T));

		//Loops while there is still space in the matrixBuffer
		while(matrixBufferPosition < matrixBufferSize){

			//Sets the starting position in the buffer
			int paddingBufferPosition = 0;

			//Loops while there is still space in the buffer
			while(paddingBufferPosition < paddingBufferSize){
				int count = 0;

				//Puts one row of data into the buffer while there is space
				//in the buffer and data left
				while(count < (int)matrix.coldim() &&
				      paddingBufferPosition < paddingBufferSize &&
				      dataOffset < (int)(matrix.coldim() * matrix.rowdim())){

					//Put entry of matrix into the padding buffer
					int row = (int) ((size_t)dataOffset / matrix.coldim());
					int col = (int) ((size_t)dataOffset % matrix.coldim());
					paddingBuffer[paddingBufferPosition] = matrix.getEntry((size_t)row, (size_t)col);

					//Increment the count for the row, paddingBuffer, and matrix
					count++;
					paddingBufferPosition++;
					dataOffset++;
					elementsPadded++;
				}

				//Padds 0's until end of padded row while there is still space in the
				//buffer
				while(count < newDimX && paddingBufferPosition < paddingBufferSize){

					//Place a zero into the paddingBuffer
					paddingBuffer[paddingBufferPosition] = 0.0;

					//Increment the count for the rwo and the paddingBuffer
					count++;
					paddingBufferPosition++;
					elementsPadded++;
				}

				if(elementsPadded >= matrixBufferSize){
					break;
				}
			}

			int transferSize;

			//Transfer the paddingBuffer to the matrixBuffer
			if((matrixBufferPosition + paddingBufferSize) <= matrixBufferSize){
				transferSize = (32 * 1024 * 1024);
			}
			//Transfer the partial paddingBuffer to the matrixBuffer
			else{
				transferSize = (int) ((size_t)(matrixBufferSize - matrixBufferPosition) * sizeof(T));
			}

			cl_int tempErrcode;
			tempErrcode = clEnqueueWriteBuffer(
				commandQue,
				matrixBuffer,
				CL_TRUE,
				((size_t)matrixBufferPosition * sizeof(T)),
				(size_t)transferSize,
				paddingBuffer,
				0,
				NULL,
				NULL);

			//Increment position in matrixBuffer by the size of the paddingBuffer
			matrixBufferPosition += paddingBufferSize;
		}

		//Deallocates buffer
		free(paddingBuffer);

		return matrixBuffer;
	}

	/**
	 * @internal
	 * Dedads an OpenCL buffer and puts it back into a BlasMatrix
	 */
	template <class Field>
	template <typename T, class Operand1>
	Operand1& OpenCLMatrixDomain<Field>::oclDepadMatrix(
		cl_mem matrixBuffer,
		int matrixBufferSize,
		int outputSize,
		int newDimX,
		Operand1& matrix) const{

		//Set starting positions
		int matrixBufferPosition = 0;
		int dataOffset = 0;

		//Allocates a 32mb buffer for depadding
		T* depaddingBuffer = (T*)malloc(32 * 1024 * 1024);

		//Calculate the size of the depadding buffer in number of elements
		const int depaddingBufferSize = (32 * 1024 * 1024 / (int)sizeof(T));

		//Loops while there are still elements in the matrixBuffer
		while(dataOffset < outputSize){
			int transferSize;

			//Transfer a full depaddingBuffer worth of elements back to the host
			if((matrixBufferPosition + depaddingBufferSize) <= matrixBufferSize){
				transferSize = (32 * 1024 * 1024);
			}
			//Transfer a partial depaddingBuffer worth of elements back to the host
			else{
				transferSize = (int) ( (size_t)(matrixBufferSize - matrixBufferPosition) * sizeof(T) );
			}

			cl_int tempErrcode;
			tempErrcode = clEnqueueReadBuffer(
				commandQue,
				matrixBuffer,
				CL_TRUE,
				((size_t)matrixBufferPosition * sizeof(T)),
				(size_t)transferSize,
				depaddingBuffer,
				0,
				NULL,
				NULL);

			//Set depaddiingBuffer start position
			int depaddingBufferPosition = 0;

			//Loops while there are still elements in the depaddingBuffer
			while(depaddingBufferPosition < depaddingBufferSize){
				int count = 0;

				//Puts one row of data into the matrix while there are elements in the
				//depaddingBuffer
				while(count < (int)matrix.coldim() &&
				      depaddingBufferPosition < depaddingBufferSize &&
				      dataOffset < outputSize){

					//Put entry of depadding buffer into the matrix
					int row = (int) ((size_t)dataOffset / matrix.coldim());
					int col = (int) ((size_t)dataOffset % matrix.coldim());
					matrix.setEntry((size_t)row, (size_t)col, depaddingBuffer[depaddingBufferPosition]);

					//Increment the count for the row, depadding buffer, and matrix
					count++;
					depaddingBufferPosition++;
					dataOffset++;
				}

				//Skip over the padding zero's
				depaddingBufferPosition += (newDimX - count);
			}

			//Increment position in matrixBuffer by the size of the paddingBuffer
			matrixBufferPosition += depaddingBufferSize;
		}

		//Deallocates buffer
		free(depaddingBuffer);

		return matrix;
	}

	/**
	 * @internal
	 * Creates an empty matrix buffer on the OpenCL device from
	 * the dimensions of the matrix
	 */
	template <class Field>
	template <typename T, class Operand1>
	cl_mem OpenCLMatrixDomain<Field>::oclCreateMatrixBuffer(
		Operand1& matrix) const{

		//Calculate dimensions after padding of matrix
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDimX = (int) (((matrix.coldim() + 15) / 16) * 16);
		int newDimY = (int) (((matrix.rowdim() + 15) / 16) * 16);

		//Allocate buffer
		cl_int tempErrcode;
		cl_mem matrixBuffer = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE,
			((size_t)newDimX * (size_t)newDimY * sizeof(T)),
			0,
			&tempErrcode);

		return matrixBuffer;
	}

	/**
	 * @internal
	 * Creates a matrix buffer on the OpenCL device from the dimensions
	 * of the passed in matrix and load the contents of the matrix into
	 * the matrix buffer
	 */
	template<class Field>
	template<typename T, class Operand1>
	cl_mem OpenCLMatrixDomain<Field>::oclCreateAndLoadMatrixBuffer(
		const Operand1 &matrix) const{

		//Calculate dimensions after padding of matrix
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDimX = (int) (((matrix.coldim() + 15) / 16) * 16);
		int newDimY = (int) (((matrix.rowdim() + 15) / 16) * 16);

		//Allocate buffer
		cl_int tempErrcode;
		cl_mem matrixBuffer = clCreateBuffer(
			context,
			CL_MEM_READ_WRITE,
			((size_t)newDimX * (size_t)newDimY * sizeof(T)),
			0,
			&tempErrcode);

		//Calculate number of elements in the matrixBuffer
		int matrixBufferSize = newDimX * newDimY;

		return oclPadMatrix<T, Operand1>(
			matrixBuffer,
			matrixBufferSize,
			newDimX,
			matrix);
	}

	/**
	 * @internal
	 * Read back the contents of the matrix buffer into the matrix and
	 * return a refence to the matrix
	 */
	template <class Field>
	template <typename T, class Operand2>
	Operand2& OpenCLMatrixDomain<Field>::oclReadMatrixBuffer(
		cl_mem matrixBuffer,
		Operand2 &matrix) const{

		//Calculate dimensions after padding of matrix
		//((A.coldim() / 16) + (A.coldim() % 16 == 0 ? 0 : 1)) * 16
		int newDimX = (int) (((matrix.coldim() + 15) / 16) * 16);
		int newDimY = (int) (((matrix.rowdim() + 15) / 16) * 16);

		//Calculate number of elements in the matrixBuffer
		int matrixBufferSize = newDimX * newDimY;

		//Calculate number of elements in the matrix
		int outputSize = (int)( matrix.coldim() * matrix.rowdim() ); //BB here we assume the size of a matrix is small, not a size_t

		return oclDepadMatrix<T, Operand2>(
			matrixBuffer,
			matrixBufferSize,
			outputSize,
			newDimX,
			matrix);
	}

} //end of namespace LinBox

#endif // __LINBOX_opencl_matrix_domain_memory_INL

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
