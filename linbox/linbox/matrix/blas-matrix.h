/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/matrix/blas-matrix.h
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __BLAS_MATRIX_H
#define __BLAS_MATRIX_H

#include <linbox/matrix/dense.h>
#include <linbox/matrix/dense-submatrix.h>

namespace LinBox {

	template <class Matrix>
	class BlasMatrix : public DenseSubMatrix<typename Matrix::Element> {


	public:
		typename Matrix::Element Element;

	private:        

		Element* _ptr; 
		size_t _stride;
		
	public:

		BlasMatrix (int m, int n) :
			DenseSubMatrix(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n), _stride(n) {}

		// Constructor from matrix (copying data)
		BlasMatrix (const Matrix& A): 	
			DenseSubMatrix(*(new DenseMatrixBase<Element> (m,n)),0,0,A.rowdim(),A.coldim()), _stride(M.coldim()) 
		{
			_ptr = _M.FullIterator();
		
			typename Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
		
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index)
				M.setEntry(iter_index.rowIndex(), iter_index.colIndex(), *iter_value);      
       
		}

		// Constructor from matrix (copying data)
		BlasMatrix (const Matrix& A, const size_t i0,const size_t j0,const size_t m, const size_t n) :
			DenseSubMatrix(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n), _stride(n) 
		{
		
			_ptr = _M.FullIterator();
			typename Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
		
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
				size_t i,j;
				i=iter_index.rowIndex();
				j=iter_index.colIndex();
				if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
					_M.setEntry(i-i0, j-j0, *iter_value);  
			}
		}
		

		// Copy Constructor of a matrix (copying data)
		BlasMatrix (const BlasMatrix<Matrix>& A) :
			DenseSubMatrix(*(new DenseMatrixBase<Element> (A._M)),0,0,A.rowdim(),A.coldim()), _stride(A._stride) {
			_ptr = _M.FullIterator();
		}

		// Copy Contructor of a matrix or submatrix (no copy is done, just through pointer)
		BlasMatrix(BlasMatrix<Matrix>& A, const size_t i=0, const size_t j=0, const size_t m=A.rowdim(), const size_t n=A.coldim()) :
			DenseSubMatrix(A,i,j,m,n), _stride(A._stride), _ptr(A._ptr+ i*A._stride+j) {}

 
		Element* getPointer() {return _ptr;}

		size_t getStride() const {return _stride;}	

	}; //end of class BlasMatrix


	// TAG for triangular blas matrix
	class BlasTag {
	public:
		typedef enum{low,up} uplo;
		typedef enum{unit,nonunit} diag;
	};


	// class of triangular blas matrix
	template <class Matrix>
	class TriangularBlasMatrix: public BlasMatrix<Matrix> {

	protected:
		
		BlasMatrix<Matrix>       &_M;
		BlasTag::uplo           _uplo;
		BlasTag::diag          _diag;
	public:

		TriangularBlasMatrix (const BlasMatrix<Matrix>& A, BlasTag::uplo x=up, BlasTag::diag y=nonunit)
			: _M(*(new BlasMatrix<Matrix>(A))) , _uplo(x), _diag(y) {}

		TriangularBlasMatrix (BlasMatrix<Matrix>& A, BlasTag::uplo x=up, BlasTag::diag y=nonunit)
			: _M(A), _uplo(x), _diag(y) {}

		BlasTag::uplo getUpLo() { return _uplo;}

		BlasTag::diag getDiag() { return _diag;}

		

	}; //end of class TriangularBlasMatrix


} //end of namespace LinBox

#endif
