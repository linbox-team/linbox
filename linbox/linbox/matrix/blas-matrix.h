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
#include <linbox/util/debug.h>


namespace LinBox {

	template <class _Element>
	class BlasMatrix : public DenseSubmatrix<_Element> {
		
	public:
		typedef _Element Element;

	protected:        
		size_t   _stride;
		bool      _alloc;
		Element    *_ptr; 		
		
	public:

		BlasMatrix () {}

		BlasMatrix (int m, int n) 
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n), _stride(n), _alloc(true) { _ptr = _M.FullIterator(); }

		// Constructor from matrix (copying data)
		template <class Matrix>
		BlasMatrix (const Matrix& A)	
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (A.coldim(),A.rowdim())),0,0,A.rowdim(),A.coldim()), _stride(A.coldim()) , _alloc(true)
		{
			_ptr = _M.FullIterator();
		
			typename Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
		
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index)
				M.setEntry(iter_index.rowIndex(), iter_index.colIndex(), *iter_value);      
       
		}

		// Constructor from matrix (copying data)
		template <class Matrix>
		BlasMatrix (const Matrix& A, const size_t i0,const size_t j0,const size_t m, const size_t n) 
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n), _stride(n) , _alloc(true) 
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
		BlasMatrix (const BlasMatrix<Element>& A) 
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (A._M)),0,0,A.rowdim(),A.coldim()), _stride(A._stride), _alloc(true) {
			_ptr = _M.FullIterator();
		}


		// Copy Contructor of a matrix (no copy is done, just through pointer)
		BlasMatrix(BlasMatrix<Element>& A) 
			: DenseSubmatrix<Element>(A), _stride(A._stride), _alloc(false), _ptr(A._ptr) {}


		// Copy Contructor of a submatrix (no copy is done, just through pointer)
		BlasMatrix(BlasMatrix<Element>& A, const size_t i, const size_t j, const size_t m, const size_t n) 
			: DenseSubmatrix<Element>(A,i,j,m,n), _stride(A._stride), _alloc(false), _ptr(A._ptr+ i*A._stride+j) {}


		
		~BlasMatrix ()  {			
			if (_alloc) {
				delete &_M;
			}
		}


		// operator = (copying data)
		BlasMatrix<Element>& operator= (const BlasMatrix<Element>& A) {		       		
			
			if (_alloc) {
				delete &_M; 
			}
			_M       = *(new DenseMatrixBase<Element> (A._M));
			_beg_row = A._beg_row;
			_end_row = A._end_row;
			_beg_col = A._beg_col;
			_end_col = A._end_col;
			_ptr     = _M.FullIterator();
			_alloc   = true;
			_stride  = A._stride;			
			return *this;
		}	
		
		Element* getPointer() const  {return _ptr;}

		Element* getWritePointer() {return _ptr;}

		size_t getStride() const {return _stride;}	

	}; //end of class BlasMatrix



	// TAG for triangular blas matrix
	class BlasTag {
	public:
		typedef enum{low,up} uplo;
		typedef enum{unit,nonunit} diag;
	};


	// class of triangular blas matrix
	template <class Element>
	class TriangularBlasMatrix: public BlasMatrix<Element> {

	protected:
		
		BlasTag::uplo           _uplo;
		BlasTag::diag           _diag;

	public:

		TriangularBlasMatrix (const size_t m, const size_t n, BlasTag::uplo x=up, BlasTag::diag y= BlasTag::nonunit)
			: BlasMatrix<Element>(m, n ) , _uplo(x), _diag(y) {}

		TriangularBlasMatrix (const BlasMatrix<Element>& A, BlasTag::uplo x=up, BlasTag::diag y= BlasTag::nonunit)
			: BlasMatrix<Element>(A) , _uplo(x), _diag(y) {}

		TriangularBlasMatrix (BlasMatrix<Element>& A, BlasTag::uplo x=up, BlasTag::diag y= BlasTag::nonunit)
			: BlasMatrix<Element>(A), _uplo(x), _diag(y) {}
		
		TriangularBlasMatrix (const TriangularBlasMatrix<Element>& A)
			: BlasMatrix<Element>(A.rowdim(),A.coldim()), _uplo(A._uplo), _diag(A._diag) {
			switch (A._uplo) {
			case BlasTag::up:
				{
					for (size_t i=0;i<A.rowdim();++i)
						for (size_t j=i;j<A.coldim();++j)
							this->setEntry(i,j,A.getEntry(i,j));
					break;
				}
			case BlasTag::low:
				{
					for (size_t i=0;i<A.rowdim();++i)
						for (size_t j=0;j<=i;++j)
							this->setEntry(i,j,A.getEntry(i,j));
					break;
				}
			default:
				throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
			}
		}		

		BlasTag::uplo getUpLo() const { return _uplo;}

		BlasTag::diag getDiag() const { return _diag;}	

	}; //end of class TriangularBlasMatrix

	template <class Element>
	struct MatrixTraits< BlasMatrix<Element> >
	{ 
		typedef BlasMatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory; 
	};

	class BlasPermutation {
		
		
	public:
		
		BlasPermutation() {};

		BlasPermutation( const size_t n ): _P(n), _order( n ) {};

		BlasPermutation( const std::vector<size_t> P ) : _P( P ), _order( P.size() ) {};

		BlasPermutation( const BlasPermutation& P): _P( P._P ), _order( P._order ) {};

		BlasPermutation& operator=( const BlasPermutation& P ){
			_P = P._P;
			_order = P._order;
			return *this;
		}
		
		
		const size_t* getPointer() const  { return &_P[0]; }
		
		size_t* getWritePointer()  { return &_P[0]; }
	
		const size_t  getOrder()  const { return _order; }

			
	
	
	protected:
		
		std::vector<size_t>  _P;
		size_t               _order;

	}; // end of class BlasPermutation

	template< class Matrix >
	class TransposedBlasMatrix {

	public:
		
		TransposedBlasMatrix ( Matrix& M ) :  _M(M) {}
		
		Matrix& getMatrix() { return _M; }
	
	protected:
		Matrix& _M;
	};
	
	template<>
	template< class Matrix >
	class TransposedBlasMatrix< TransposedBlasMatrix< Matrix > > : public Matrix {
		
	public:
		TransposedBlasMatrix ( Matrix& M ) :  Matrix(M){}	
		TransposedBlasMatrix ( const Matrix& M ) :  Matrix(M){}	
		
	};
	
} //end of namespace LinBox

#endif
