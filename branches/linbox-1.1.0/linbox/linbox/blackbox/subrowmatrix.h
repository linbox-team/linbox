/* -*- mode: C++ -*- */
/* Author: Zhendong wan
*/

#ifndef __SUBROWMATRIX_H__
#define __SUBROWMATRIX_H__

#include <iostream>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/blackbox/blackbox-interface.h>

namespace LinBox {
	/** 
	 * \brief submatrix consisting contiguous rows of a row based matrix.

\ingroup blackbox
	 * submatrix consisting of rows [i1, ..., i2] of a row based matrix.
	 */
	//@{
	template<class Matrix,
	   	 class MatrixCategory = typename MatrixTraits<Matrix>::MatrixCategory>

	class SubRowMatrix: public BlackboxInterface ;

	/** matrix with RowMatrixTag
	  * Support Row, apply, applyTranspose
	  */
	template <class Matrix>
	class SubRowMatrix<Matrix, MatrixCategories::RowMatrixTag> 
	: public BlackboxInterface {

	public:
		typedef typename Matrix::Field Field;
		  
	protected:

		/* try to use const Matrix* 
		  * but getting too much compiling error messages.
		  * Maybe do it in future.
		  * Z. Wan
		  */
		Matrix* _BB;
		
		size_t _row;
		
		size_t _rowdim;

		MatrixDomain<Field> _MD;

	public:

	 	typedef typename MatrixCategories::RowMatrixTag MatrixCategory;

		typedef typename Field::Element Element;

		typedef typename Matrix::RowIterator            RowIterator;
		typedef typename Matrix::ConstRowIterator       ConstRowIterator;
		typedef typename Matrix::Row                    Row;
		typedef typename Matrix::ConstRow               ConstRow;

		typedef typename Matrix::ConstRawIterator RawIterator;
		typedef typename Matrix::ConstRawIterator ConstRawIterator;
		typedef typename Matrix::ConstRawIndexedIterator RawIndexedIterator;
		typedef typename Matrix::ConstRawIndexedIterator ConstRawIndexedIterator;

		SubRowMatrix (const Matrix* BB,
			   size_t row,
			   size_t rowdim) 
			: _row(row), _rowdim(rowdim), _MD(BB->field()) {

			_BB = const_cast<Matrix*&>(BB);

			}


		virtual ~SubRowMatrix ( ) {}

		size_t rowdim ( ) const {

			return _rowdim;
		}

		size_t coldim ( ) const {

			return _BB -> coldim();
		}

		const Field& field( ) const {

			return _BB -> field ();
		}
		
		template <class OutVector, class InVector>
		OutVector &apply (OutVector &y, const InVector &x) const {
			
			
			return _MD.vectorMul (y, *this, x); 
		}
		
		
		template <class OutVector, class InVector>
		OutVector &applyTranspose (OutVector& y, const InVector &x) const {

			return _MD.vectorMul (y, TransposeMatrix<SubRowMatrix> (*this), x);
		}



            template<typename _Tp1> 
            struct rebind                           
            { typedef SubRowMatrix<Matrix::rebind<_Tp1>::other, MatrixCategory> other; };



		RowIterator rowBegin () {

			return _BB -> rowBegin() + _row;

		}
		
		RowIterator rowEnd () {

			return _BB -> rowBegin() + (_row + rowdim);
		}

		
		ConstRowIterator rowBegin () const {
			
			return _BB -> rowBegin() + _row;
		}
			
		ConstRowIterator rowEnd () const {
			
			return _BB -> rowBegin() + (_row + _rowdim);
		}


		std::ostream& write (std::ostream& out) const {

#ifdef DEBUG

			std::cout << "Row, Col " << rowdim () << ' ' << coldim () << '\n';
#endif

			ConstRowIterator row_p;

			VectorDomain<Field> vd (field());


			for (row_p = rowBegin(); row_p != rowEnd(); ++ row_p) {

				vd. write (out, *row_p);
			
				std::cout << '\n';

			}

			return out;

		}

		/* I do not need them currently now, 
		  * I maybe implement in future
		  */
		//RawIterator rawBegin ();
		//RawIterator rawEnd ();
		//ConstRawIterator rawBegin () const;
		//ConstRawIterator rawEnd () const;
	};
	//@}
}
#endif
