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

#ifndef __LINBOX_blas_matrix_H
#define __LINBOX_blas_matrix_H


#include <linbox/solutions/getentry.h>
#include <linbox/matrix/dense.h>
#include <linbox/matrix/dense-submatrix.h>
#include <linbox/util/debug.h>
#include <linbox/matrix/matrix-category.h>

namespace LinBox 
{

	template<class Element>
	class BlasMatrix;
	
	template <class Element>
	class MatrixContainerTrait<BlasMatrix<Element> > 
	{
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};
	
	template <class Element>
	class MatrixContainerTrait<const BlasMatrix<Element> > 
	{
	public:
		typedef MatrixContainerCategory::BlasContainer Type;
	};

	// @brief Limited docs so far.
	template <class _Element>
	class BlasMatrix : public DenseSubmatrix<_Element> 
	{
		
	public:
            typedef _Element Element;
            typedef typename DenseSubmatrix<_Element>::Element * pointer;
            typedef BlasMatrix<Element> Self_t;

	protected:        
		size_t   _stride;
		bool      _alloc;
		pointer  _ptr; 		
		 
	public:
		typedef typename DenseSubmatrix<_Element>::RawIterator RawIterator ; // for nvcc

		BlasMatrix ()
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (0,0)),0,0,0,0), _stride(0),  _alloc(true) { _ptr = this->_M->FullIterator(); }

		
		
		BlasMatrix (int m, int n) 
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n), _stride(n), _alloc(true) { _ptr = this->_M->FullIterator();}

		BlasMatrix (size_t m, size_t n) 
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (m,n)),0,0,m,n), _stride(n), _alloc(true) { _ptr = this->_M->FullIterator();}



		/** Constructor from a matrix stream */
		template< class Field >
		BlasMatrix(MatrixStream<Field>& ms)
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (ms))),  _alloc(true) { _ptr = this->_M->FullIterator(); _stride= this->coldim();}


		// Generic copy constructor from either a blackbox or a matrix container
		template <class Matrix>
		BlasMatrix (const Matrix &A)
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (A.rowdim(),A.coldim())),0,0,A.rowdim(),A.coldim()), _stride(A.coldim()) , _alloc(true)
		{
			_ptr = this->_M->FullIterator();
			createBlasMatrix(A, typename MatrixContainerTrait<Matrix>::Type());
		}
		
		// Generic copy constructor from either a blackbox or a matrix container (allow submatrix)
		template <class Matrix>
		BlasMatrix (const Matrix& A, const size_t i0,const size_t j0,const size_t m, const size_t n)
			: DenseSubmatrix<Element>( *(new DenseMatrixBase<Element> (A.rowdim(),A.coldim())),0,0,A.rowdim(),A.coldim()), _stride(A.coldim()) , _alloc(true)
		{
			_ptr = this->_M->FullIterator();
			createBlasMatrix(A, i0, j0, m, n, typename MatrixContainerTrait<Matrix>::Type());
		}

#if 0
		template<class _Matrix, class _Field>
		BlasMatrix (const _Matrix &A,  const _Field &F) 
			: DenseSubmatrix<Element>( *(new DenseMatrixBase<Element> (A.rowdim(),A.coldim())),0,0,A.rowdim(),A.coldim() ), 
			_stride(A.coldim()) , 
			_alloc(true)
		{
			_ptr = this->_M->FullIterator() ;
			rebind<_Field>()(*this,A,F);

		}
#endif

		
		// Copy data according to blas container structure
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A, MatrixContainerCategory::BlasContainer)	
			
		{
			typename _Matrix::ConstRawIterator         iter_value = A.rawBegin();
			RawIterator  iter_addr = this->rawBegin();			
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_addr)
				*iter_addr = *iter_value;			
		}

		// Copy data according to Matrix container structure
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Container)	
			
		{
			// With both iterators, it is Segfaulting !!!!		
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
			for (;iter_index != A.rawIndexedEnd(); ++iter_index)
				this->_M->setEntry( iter_index.rowIndex(), 
                                                    iter_index.colIndex(), 
                                                    A.getEntry(iter_index.rowIndex(),iter_index.colIndex())
                                                    );
		}
		
		
		// Copy data according to blackbox structure
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Blackbox)	
			
		{
			typedef typename Matrix::Field Field;
			typename Field::Element one, zero;
			Field F = A.field();
			F. init(one, 1);
			F. init(zero, 0);
			
			std::vector<typename Field::Element> e(A.coldim(), zero), tmp(A.rowdim());
			typedef typename DenseSubmatrix<_Element>::ColIterator ColIterator ;
			ColIterator col_p;
			
			typename BlasMatrix<Element>::Col::iterator elt_p;
			
			typename std::vector<typename Field::Element>::iterator e_p, tmp_p;
			
			
			//for (col_p = colBegin(), e_p = e.begin();
			for (col_p = DenseSubmatrix<_Element>:: colBegin(), e_p = e.begin();
			     e_p != e.end(); ++ col_p, ++ e_p) {
				
				F.assign(*e_p, one);
				
				A.apply (tmp, e);
				
				for (tmp_p = tmp.begin(), elt_p = col_p -> begin();
				     tmp_p != tmp.end(); ++ tmp_p, ++ elt_p)
					
					F.assign(*elt_p, *tmp_p);
				
				F.assign(*e_p, zero);
			}			
		}

		// Copy data according to Matrix container structure (allow submatrix)
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A, const size_t i0,const size_t j0,const size_t m, const size_t n, MatrixContainerCategory::Container)
		{
		
			typename _Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename _Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
		
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
				size_t i,j;
				i=iter_index.rowIndex();
				j=iter_index.colIndex();
				if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
					this->_M->setEntry(i-i0, j-j0, *iter_value);  
			}
		}

		// Copy data according to Matrix container structure (allow submatrix)
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, const size_t i0,const size_t j0,const size_t m, const size_t n, MatrixContainerCategory::BlasContainer)
		{
		
			typename Matrix::ConstRawIterator         iter_value = A.rawBegin();
			typename Matrix::ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
		
			for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
				size_t i,j;
				i=iter_index.rowIndex();
				j=iter_index.colIndex();
				if (( i >= i0) && (i< i0+m) && (j >= j0) && (j < j0+n))
					this->_M->setEntry(i-i0, j-j0, *iter_value);  
			}
		}


		// Copy data according to blackbox structure (allow submatrix)
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, const size_t i0,const size_t j0,const size_t m, const size_t n, MatrixContainerCategory::Blackbox) 			
		{
			std::cerr << __func__ << ": not implemented yet" << std::flush << std::endl; 
			exit(-1) ;
			//! @todo need to be implemented by succesive apply
			//! @todo lancer une exception générique "not implemented yet"
		}


#if 0
 		// Constructor from matrix (no copy)
 		BlasMatrix ( DenseMatrixBase<Element>& A )	
 			: DenseSubmatrix<Element>(A,0,0,A.rowdim(),A.coldim()), _stride(A.coldim()) , _alloc(false)
 		{ _ptr = this->_M->FullIterator();}

 		// Constructor from matrix (no copy )
 		BlasMatrix ( DenseMatrixBase<Element>& A, const size_t i0,const size_t j0,const size_t m, const size_t n) 
 			: DenseSubmatrix<Element>(A,i0,j0,m,n), _stride(A.coldim()) , _alloc(false) 
 		{_ptr = this->_M->FullIterator();}
#endif


		// Copy Constructor of a matrix (copying data)
		BlasMatrix (const BlasMatrix<Element>& A) 
			: DenseSubmatrix<Element>(*(new DenseMatrixBase<Element> (*A._M)),0,0,A.rowdim(),A.coldim()), _stride(A._stride), _alloc(true) 
		{ _ptr = this->_M->FullIterator(); }
		
#if 0
 		// Copy Contructor of a matrix (no copy is done, just through pointer)
 		BlasMatrix(BlasMatrix<Element>& A) 
 			: DenseSubmatrix<Element>(A), _stride(A._stride), _alloc(false), _ptr(A._ptr) {}


 		// Copy Contructor of a submatrix (no copy is done, just through pointer)
 		BlasMatrix(BlasMatrix<Element>& A, const size_t i, const size_t j, const size_t m, const size_t n) 
 			: DenseSubmatrix<Element>(A,i,j,m,n), _stride(A._stride), _alloc(false), _ptr(A._ptr+ i*A._stride+j) {}
#endif

		
		~BlasMatrix ()  
		{			
			if (_alloc) 
				delete this->_M;
		}

                // Rebind operator
		template<typename _Tp1>
		struct rebind
		{ 
			typedef BlasMatrix<typename _Tp1::Element> other; 

			void operator() (other & Ap, const Self_t& A, const _Tp1& F) 
			{

				typedef typename DenseSubmatrix<_Element>::ConstRawIndexedIterator ConstRawIndexedIterator ;
				typedef typename DenseSubmatrix<_Element>::ConstRawIterator ConstRawIterator ;
				ConstRawIterator         iter_value = A.rawBegin();
				ConstRawIndexedIterator  iter_index = A.rawIndexedBegin();
				typename _Tp1::Element tmp;
				for (;iter_value != A.rawEnd(); ++iter_value,++iter_index){
					F.init(  tmp, *iter_value ); 
					Ap.setEntry(iter_index.rowIndex(), iter_index.colIndex(),tmp);
				}
			}
		};


		// operator = (copying data)
		BlasMatrix<Element>& operator= (const BlasMatrix<Element>& A) 
		{	

			DenseMatrixBase<Element> *tmp= this->_M;
			this->_M       = new DenseMatrixBase<Element>(*A._M);
			if (_alloc) {
				delete tmp; 
			}
			this->_beg_row = A._beg_row;
			this->_end_row = A._end_row;
			this->_beg_col = A._beg_col;
			this->_end_col = A._end_col;
			_ptr     = this->_M->FullIterator();
			_alloc   = true;
			_stride  = A._stride;			

			return *this;
		}

		pointer getPointer() const  {return _ptr;}
		pointer& getWritePointer() {return _ptr;}

		size_t getStride() const {return _stride;}	
		size_t& getWriteStride() {return _stride;}	

	}; // end of class BlasMatrix


	// TAG for triangular blas matrix
	class BlasTag 
	{
	public:
		typedef enum{low,up} uplo;
		typedef enum{unit,nonunit} diag;
	};


	// class of triangular blas matrix
	template <class Element>
	class TriangularBlasMatrix: public BlasMatrix<Element> 
	{

	protected:
		
		BlasTag::uplo           _uplo;
		BlasTag::diag           _diag;

	public:

		TriangularBlasMatrix (const size_t m, const size_t n, BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit)
			: BlasMatrix<Element>(m, n ) , _uplo(x), _diag(y) {}

		TriangularBlasMatrix (const BlasMatrix<Element>& A, BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit)
			: BlasMatrix<Element>(A) , _uplo(x), _diag(y) {}

		TriangularBlasMatrix (BlasMatrix<Element>& A, BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit)
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
                                    for (size_t i=0;i<A.rowdim();++i) {
                                        for (size_t j=0;j<=i;++j)
                                            this->setEntry(i,j,A.getEntry(i,j));
                                    }
                                    
                                    break;
				}
			default:
				throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
			}
		}		

            template<class Matrix>
            TriangularBlasMatrix (const Matrix& A, BlasTag::uplo x=BlasTag::up, BlasTag::diag y= BlasTag::nonunit)
                    : BlasMatrix<Element>(A.rowdim(),A.coldim()), _uplo(x), _diag(y) {
                switch (x) {
                    case BlasTag::up:
                    {
                        for (size_t i=0;i<A.rowdim();++i){
                            for (size_t j=i;j<A.coldim();++j) {
                                Element tmp;
                                this->setEntry(i,j,getEntry(tmp, A,i,j));
                            }
                        }
                        break;
                    }
                    case BlasTag::low:
                    {
                        for (size_t i=0;i<A.rowdim();++i) {
                            for (size_t j=0;j<=i;++j) {
                                Element tmp;
                                this->setEntry(i,j,getEntry(tmp, A,i,j));
                            }
                        }
                        
                        break;
                    }
                    default:
                        throw LinboxError ("Error in copy constructor of TriangularBlasMatrix (incorrect argument)");
                }
            }		

		BlasTag::uplo getUpLo() const { return _uplo;}

		BlasTag::diag getDiag() const { return _diag;}	

	}; // end of class TriangularBlasMatrix

	template <class Element>
	struct MatrixTraits< BlasMatrix<Element> >
	{ 
		typedef BlasMatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory; 
	};

	template <class Element>
	struct MatrixTraits< const BlasMatrix<Element> >
	{ 
		typedef const BlasMatrix<Element> MatrixType;
		typedef typename MatrixCategories::RowColMatrixTag MatrixCategory; 
	};


	/** Class used for permuting indices. For example, create a vector (0 1 2 ...) over size_t,
	 *  then apply a permutation to it using a BlasMatrixDomain to get the natural representation of the permutation.
	 */
	class indexDomain 
	{
	public:
		typedef size_t Element;
	public:
		indexDomain() {};
		template <class ANY>
		size_t init(size_t& dst, const ANY& src) const {
			return dst = static_cast<size_t>(src);
		}
		template <class ANY>
		size_t assign(ANY& dst, const size_t& src) const {
			return dst = static_cast<ANY>(src);
		}
	};


// Dan Roche 7-8-04 Changed _P to _PP to avoid confict with a macro defined in
// <iostream> somewhere.
	class BlasPermutation 
	{
		
		
	public:
		
		BlasPermutation() {};

		BlasPermutation( const size_t n ): _PP(n), _order( n ) {};

		BlasPermutation( const std::vector<size_t> P ) : _PP( P ), _order( P.size() ) {};

		BlasPermutation( const BlasPermutation& P): _PP( P._PP ), _order( P._order ) {};

		BlasPermutation& operator=( const BlasPermutation& P ){
			_PP = P._PP;
			_order = P._order;
			return *this;
		}
		
		
		const size_t* getPointer() const  { return &_PP[0]; }
		
		size_t* getWritePointer()  { return &_PP[0]; }
	
		//const size_t  getOrder()  const { return _order; } // BB: "warning: type qualifier on return type is meaningless"
		size_t  getOrder()  const { return _order; }

		BlasPermutation& extendTrivially(const size_t newSize) {
			if (newSize < _order) 
				std::cerr << "WARNING: attempting to reduce size of permutation.";
			_PP.resize(newSize);
			for (size_t i=_order; i<newSize; i++)
				_PP[i] = i;
			_order = newSize;
			return *this;
		};
	
	protected:
		
		std::vector<size_t>  _PP;
		size_t               _order;

	}; // end of class BlasPermutation

	template< class Matrix >
	class TransposedBlasMatrix 
	{

	public:
		
		TransposedBlasMatrix ( Matrix& M ) :  _M(M) {}
		
		Matrix& getMatrix() const { return _M; }
	
	protected:
		Matrix& _M;
	};
	
#if !defined(__INTEL_COMPILER) && !defined(__CUDACC__)
	template <>
#endif
	template< class Matrix >
	class TransposedBlasMatrix< TransposedBlasMatrix< Matrix > > : public Matrix 
	{
		
	public:
		TransposedBlasMatrix ( Matrix& M ) :  Matrix(M){}	
		TransposedBlasMatrix ( const Matrix& M ) :  Matrix(M){}	
		
	};
	

} // end of namespace LinBox

#endif // __LINBOX_blas_matrix_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
