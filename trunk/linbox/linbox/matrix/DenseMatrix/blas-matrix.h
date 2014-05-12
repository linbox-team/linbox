/*
 * Copyright (C) 2004 Pascal Giorgi, Clément Pernet
 *               2013, 2014 the LinBox group
 *
 * Written by :
 *               Pascal Giorgi  pascal.giorgi@ens-lyon.fr
 *               Clément Pernet clement.pernet@imag.fr
 *               Brice Boyer    <bboyer@imag.fr>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file matrix/DenseMatrix/blas-matrix.h
 * @ingroup densematrix
 * A \c BlasMatrix<\c _Field > represents a matrix as an array of
 * <code>_Field::Element</code>s. It also has the BlasBlackbox interface.
 *
 */

#ifndef __LINBOX_matrix_densematrix_blas_matrix_H
#define __LINBOX_matrix_densematrix_blas_matrix_H

#include "linbox/util/debug.h"
#include "linbox/linbox-tags.h"
#include "linbox/vector/stream.h"
#include "linbox/field/hom.h"
#include "linbox/vector/blas-vector.h"

#include "linbox/matrix/matrix-category.h"
#include "linbox/matrix/matrix-traits.h"
#include "linbox/util/matrix-stream.h"

#include "linbox/field/modular.h" // just for checkBlasApply
#include "linbox/field/modular-balanced.h" // just for checkBlasApply

//! @bug this does not belong here.
#include "blas-transposed-matrix.h"
#include "linbox/matrix/MatrixDomain/matrix-domain.h"
#include "linbox/matrix/MatrixDomain/apply-domain.h"

namespace LinBox
{ /*  not generic wrt Field (eg NTL_ZZ_p) */
	namespace Protected
	{

		//!@bug this does not seem right for float or any non M/modular field: doing blas wherever we have a fflas-ffpack field (?)
		//! @bug should return true for some UnparametricField
		template <class Field>
		bool checkBlasApply(const Field &F, size_t n)
		{

			return false;
			// integer chara, card;
			// F.characteristic(chara);
			// F.cardinality(card);

			// if ((chara != card) || chara == 0)
				// return false;
			// else
				// if (n*chara*chara < integer("9007199254740992"))
					// return true;
				// else
					// return false;
		}

		template<>
		bool checkBlasApply(const Modular<double> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<double> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<float> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<float> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<double>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<double>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<float>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<float>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<int64_t> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<int64_t> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<int32_t> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<int32_t> &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<int64_t>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<int64_t>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const Modular<int32_t>::Father_t &, size_t)
		{
			return true;
		}

		template<>
		bool checkBlasApply(const ModularBalanced<int32_t>::Father_t &, size_t)
		{
			return true;
		}

	}
}

namespace LinBox
{ /*  Blas Matrix */
	template<class Matrix>
	class MatrixDomain;

	/*! Dense matrix representation.
	 * @ingroup matrix
	 * A \p BlasMatrix is a matrix of \p _Field::Element, with the structure of BLAS matrices.
	 * It is basically a vector of \p _Field::Element.
	 * In the Mother model, a \p BlasMatrix is allocated by the user.
	 *@bug why not BlasMatrixDomain ?
	 */
	template <class _Field, class _Storage>
	class BlasMatrix {
		// private :

	public:
		typedef _Field                                  Field;
		typedef typename Field::Element               Element;    //!< Element type
		typedef _Storage                                  Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef typename Rep::pointer                 pointer;    //!< pointer type to elements
		typedef const pointer                   const_pointer;    //!< const pointer type
		typedef BlasMatrix<Field,Rep>                  Self_t;    //!< Self typeype
		typedef const BlasMatrix<Field,Rep>       constSelf_t;    //!< Self typeype

                typedef BlasSubmatrix<Self_t>           subMatrixType;    //!< Submatrix type
		typedef BlasSubmatrix<constSelf_t> constSubMatrixType;    //!< Submatrix type
                typedef Self_t                             matrixType;    //!< matrix type
                typedef constSelf_t                   constMatrixType;    //!< matrix type
                typedef Self_t                               blasType;    //!< blas matrix type

	protected:
		size_t			    _row;
		size_t			    _col;
		Rep			    _rep;
	public:
		bool		     _use_fflas ; //! @bug why public ?
	protected:
		pointer			    _ptr;
	public:
	// protected:
	       	const Field		    * _field; //! @bug why public ?
		MatrixDomain<Field>    _MD; //! @bug why public ?
		VectorDomain<Field>    _VD;
		// applyDomain<subMatrixType>    _AD; //! @bug why public ?
		// applyDomain<Self_t>    _AD; //! @bug why public ?


	private:

#if 0
		void makePointer()
		{
#if 0
			if (_row && _col) {
				_ptr = malloc( _row*_col*sizeof(_Element) ) ;
				linbox_check(_ptr);
			}
			else
				_ptr = NULL ;
#endif
			_rep = Rep(_row*_col);
			_ptr = &_rep[0];
		}
#endif

		/*! @internal
		 * @name create BlasMatrix
		 * @{ */

		/*! @internal
		 * Copy data according to blas container structure.
		 * Specialisation for BlasContainer.
		 */
		void createBlasMatrix (const Self_t & A) ;

		/*! @internal
		 * Copy data according to blas container structure.
		 * Specialisation for BlasContainer.
		 */
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A, MatrixContainerCategory::BlasContainer) ;

		/*! @internal
		 * Copy data according to Matrix container structure.
		 * Specialisation for Container
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Container) ;

		/*! @internal
		 * Copy data according to blackbox structure.
		 * Specialisation for Blackbox.
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A, MatrixContainerCategory::Blackbox) ;

		/*! @internal
		 * Copy data according to Matrix container structure (allow submatrix).
		 * Specialisation for Container
		 */
		template <class _Matrix>
		void createBlasMatrix (const _Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::Container) ;

		/*! @internal
		 * Copy data according to Matrix container structure (allow submatrix).
		 * Specialisation for BlasContainer.
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::BlasContainer) ;

		/*! @internal
		 * Copy data according to blackbox structure (allow submatrix).
		 * Specialisation for Blackbox matrices
		 * @todo need to be implemented by succesive apply
		 */
		template <class Matrix>
		void createBlasMatrix (const Matrix& A,
				       const size_t i0,const size_t j0,
				       const size_t m, const size_t n,
				       MatrixContainerCategory::Blackbox) ;

		/*!@internal constructor from vector of elements.
		 * @param v pointer to \c Element s
		 */
		void createBlasMatrix ( const Element * v) ;

		/*!@internal constructor from vector of elements.
		 * @param v std::vector of \c Element s
		 */
		void createBlasMatrix ( const std::vector<Element> & v) ;
		/*! @internal
		 * @}
		 */

	public:

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*! Allocates a new \f$ 0 \times 0\f$ matrix (shaped and ready).
		*/
		BlasMatrix (const _Field &F) ;

		// /*! Allocates a new bare \f$ 0 \times 0\f$ matrix (unshaped, unready).
		// */
		// BlasMatrix () ;

		/// (Re)allocates a new \f$ m \times n\f$ zero matrix (shaped and ready).
		void init(const _Field & F, const size_t & r = 0, const size_t & c = 0);

		/*! Allocates a new \f$ m \times n\f$ zero matrix (shaped and ready).
		 * @param F
		 * @param m rows
		 * @param n cols
		 */
		//@{

		BlasMatrix (const _Field &F, const size_t & m, const size_t &n) ;

		//@}


		/*! Constructor from a matrix stream.
		 * @param ms matrix stream.
		 */
		BlasMatrix(MatrixStream<_Field>& ms) ;

		/*! Generic copy constructor from either a blackbox or a matrix container.
		 * @param A matrix to be copied
		 */
		template <class Matrix>
		BlasMatrix (const Matrix &A) ;

		/*! Generic copy constructor from either a blackbox or a matrix container (allow submatrix).
		 * @param A matrix to be copied
		 * @param i0
		 * @param j0
		 * @param m rows
		 * @param n columns
		 */
		template <class Matrix>
		BlasMatrix (const Matrix& A,
			    const size_t & i0, const size_t & j0,
			    const size_t & m,  const size_t & n) ;

		/*! Constructor.
		 * @param A matrix to be copied
		 * @param F Field
		 */
		template<class _Matrix>
		BlasMatrix (const _Matrix &A,  const _Field &F) ;

		/*! Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		BlasMatrix (const Self_t & A) ;

		/*- Copy Constructor of a matrix (copying data).
		 * @param A matrix to be copied.
		 */
		// BlasMatrix (const BlasSubmatrix<Field,Rep>& A) ;

		/*! Create a BlasMatrix from a vector of elements
		 * @param F
		 * @param v
		 * @param m
		 * @param n
		 */
		BlasMatrix (const _Field &F, const std::vector<Element>& v,
			    const size_t &m , const size_t &n) ;

		/*! Create a BlasMatrix from an array of elements
		 * @param F
		 * @param v
		 * @param m
		 * @param n
		 */
		BlasMatrix (const _Field &F, const Element * v,
			    const size_t & m, const size_t & n) ;


		/** Constructor using a finite vector stream (stream of the rows).
		 * @param  F The field of entries; passed so that arithmetic may be done
		 *           on elements
		 * @param  stream A vector stream to use as a source of vectors for this
		 *                matrix
		 */
		template <class StreamVector>
		BlasMatrix (const Field &F, VectorStream<StreamVector> &stream) ;

		/// Destructor.
		~BlasMatrix () ;

		//! operator = (copying data)
		Self_t& operator= (const Self_t& A) ;

		//! make sure we actually copy
		template<class Matrix>
		BlasMatrix &copy( const Matrix & B)
		{
			for (size_t i = 0 ; i < rowdim() ; ++i)
				for (size_t j = 0 ; j < coldim() ; ++j) {
					setEntry(i,j,B.getEntry(i,j));
				}
			return *this;

		}


		//! Rebind operator
		template<typename _Tp1>
		struct rebind ;

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the number of rows in the matrix.
		 * @returns Number of rows in matrix
		 */
		size_t rowdim() const ;

		/** Get the number of columns in the matrix.
		 * @returns Number of columns in matrix
		 */
		size_t coldim() const ;

		/*! Get the stride of the matrix.
		 */
		size_t getStride() const;
		size_t stride() const { return getStride() ;}

		/*!Get a reference to the stride of the matrix.
		 * Modify stride this way.
		 */
		size_t& getWriteStride();


		/** Resize the matrix to the given dimensions.
		 * The state of the matrix's entries after a call to this method is
		 * undefined
		 * @param m Number of rows
		 * @param n Number of columns
		 * @param val
		 */
		void resize (const size_t &m, const size_t &n, const Element& val = Element()) ;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		/*! @internal
		 * Get read-only pointer to the matrix data.
		 */
		pointer getPointer() const ;

		const_pointer &getConstPointer() const ;

		Rep & refRep() { return _rep ;}
		const Rep & getRep() const { return _rep ;}


		/*! @internal
		 * Get write pointer to the matrix data.
		 * Data may be changed this way.
		 */
		pointer& getWritePointer() ;

		/** Set the entry at the (i, j) position to a_ij.
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, const Element &a_ij) ;

		/** Get a writeable reference to the entry in the (i, j) position.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @returns Reference to matrix entry
		 */
		Element &refEntry (size_t i, size_t j) ;

		/** Get a read-only reference to the entry in the (i, j) position.
		 * @param i Row index
		 * @param j Column index
		 * @returns Const reference to matrix entry
		 */
		const Element &getEntry (size_t i, size_t j) const ;

		/** Copy the (i, j) entry into x, and return a reference to x.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @returns Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const ;

		///////////////////
		// TRANSPOSE &AL //
		///////////////////

		/*! Creates a transposed matrix of \c *this.
		 * @param[in] tM
		 * @return the transposed matrix of this.
		 */
		Self_t transpose(Self_t & tM) const ;


		/*! Transpose (inplace).
		 * If rows and columns agree, we can transpose inplace.
		 */
		template<bool _IP>
		void transpose() ;

		void transpose() ;

		/*! Reverse the rows of a matrix.
		 * This is done inplace.
		 * Let J=antiDiag(1) (or the matrix of the reverse
		 * permutation or the matrix (i,j) = (i+j+1==m)). Then,
		 * we compute A <- J.A;
		 */
		void reverseRows() ;

		/*! Reverse the columns of a matrix.
		 * This is done inplace.
		 * This is A <- J.A
		 */
		void reverseCols() ;

		/*! Reverse the rows/columns of a matrix.
		 * This is done inplace.
		 * This is A <- J.A.J
		 */
		void reverse() ;

		// init to field zero elements
		void zero() ;
		// init to random field elements
		void random()
		{
			subMatrixType B(*this, 0, 0, rowdim(), coldim());
			B.random();
		}

		template<class Rand>
		void random(const Rand&)
		{
			return random();
		}

		///////////////////
		//      I/O      //
		///////////////////

		/** Read the matrix from an input stream.
		 * The stream is in SMS, DENSE, or MatrixMarket format and is autodetected.
		 * @param file Input stream from which to read
		 */
		std::istream &read (std::istream &file);

		/// Write the matrix in MatrixMarket format.
		std::ostream &write (std::ostream &os) const
		{
			// std::cout << "writing" << std::endl;
			constSubMatrixType B(*this, 0, 0, rowdim(), coldim());
			// std::cout << "......." << std::endl;
			return B.write(os);
		}

		/** Write the matrix to an output stream.
		 * @param os Output stream to which to write
		 * @param f write in some format (@ref Tag::FileFormat::Format). Default is Maple's.
		 */
		std::ostream &write (std::ostream &os,
				     LINBOX_enum (Tag::FileFormat) f/* = Tag::FileFormat::Maple*/) const
		{
			constSubMatrixType B(*this, 0, 0, rowdim(), coldim());
			return B.write(os, f);
		}

		/*! @deprecated Only for compatiblity.
		 */
		std::ostream &write (std::ostream &os,
				     bool mapleFormat) const
		{
			constSubMatrixType B(*this, 0, 0, rowdim(), coldim());
			return B.write(os, mapleFormat);
		}




		///////////////////
		//   ITERATORS   //
		///////////////////

		/** @name Column of rows iterator
		 * \brief
		 * The column of rows iterator traverses the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in dense format
		 */
		//@{
		typedef Subvector<typename Rep::iterator, typename Rep::const_iterator> Row;
		typedef Subvector<typename Rep::const_iterator>                    ConstRow;

		/*!  Row Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class RowIterator;
		/*! Const Row Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ConstRowIterator;

		RowIterator      rowBegin ();
		RowIterator      rowEnd ();
		ConstRowIterator rowBegin () const;
		ConstRowIterator rowEnd   () const;
		//@}

		/** @name Row of columns iterator
		 * \brief
		 * The row of columns iterator traverses the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 */
		//@{
		typedef Subvector<Subiterator<typename Rep::iterator> >            Col;
		typedef Subvector<Subiterator<typename Rep::const_iterator> > ConstCol;
		typedef Col           Column;
		typedef ConstCol ConstColumn;

		/*! Col Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ColIterator;
		/*! Const Col Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ConstColIterator;

		ColIterator      colBegin ();
		ColIterator      colEnd ();
		ConstColIterator colBegin () const;
		ConstColIterator colEnd ()   const;
		//@}

		/** @name Iterator
		 * \brief
		 *
		 * The iterator is a method for accessing all entries in the matrix
		 * in some unspecified order. This can be used, e.g. to reduce all
		 * matrix entries modulo a prime before passing the matrix into an
		 * algorithm.
		 */
		//@{
		typedef typename Rep::iterator Iterator;
		typedef typename Rep::const_iterator ConstIterator;

		Iterator      Begin ();
		Iterator      End   ();
		ConstIterator Begin () const;
		ConstIterator End   () const;
		//@}

		/** @name Raw Indexed iterator
		 * \brief
		 *
		 * Like the raw iterator, the indexed iterator is a method for
		 * accessing all entries in the matrix in some unspecified order.
		 * At each position of the the indexed iterator, it also provides
		 * the row and column indices of the currently referenced entry.
		 * This is provided through it's \c rowIndex() and \c colIndex() functions.
		 */
		//@{
		class IndexedIterator;
		/*! Const Indexed Iterator.
		 * @ingroup iterators
		 * @brief NO DOC
		 */
		class ConstIndexedIterator;

		IndexedIterator      IndexedBegin ();
		IndexedIterator      IndexedEnd   ();
		ConstIndexedIterator IndexedBegin () const;
		ConstIndexedIterator IndexedEnd   () const;
		//@}

		/** Retrieve a reference to a row.
		 * Since rows may also be indexed, this allows A[i][j] notation
		 * to be used.
		 * @param i Row index
		 * @bug Rows and Cols should be BlasVectors
		 */
		//@{
		Row      operator[] (size_t i) ;
		ConstRow operator[] (size_t i) const ;
		//@}

		///////////////////
		//     MISC     //
		///////////////////


		/** Compute column density.
		 * @param v
		 */
		template <class Vector>
		Vector &columnDensity (Vector &v) const ;

		size_t size()
		{
			return _row * _col;
		}

                void finalize() {}

		///////////////////
		//   BLACK BOX   //
		///////////////////


		template <class Vector1, class Vector2>
		Vector1&  apply (Vector1& y, const Vector2& x) const ;

		template<class _VRep>
		BlasVector<Field,_VRep>&  apply (BlasVector<Field,_VRep>& y, const BlasVector<Field,_VRep>& x) const ;

		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const ;

		const _Field& field() const;
		//_Field& field() ;
		// void setField(const _Field & F) { _field = F ; };

		template<class uselessTag>
		void changeFieldSpecialised( _Field & G,
					     // MatrixDomain<_Field> & MD,
					     VectorDomain<_Field> & VD,
					     const _Field & F,
					     const uselessTag & m)
		{
			// don't do anything (?)
			return;
		}

		void changeFieldSpecialised(      _Field & G,
						  // MatrixDomain<_Field> & MD,
						  VectorDomain<_Field> & VD,
					    const _Field & F,
					    const RingCategories::ModularTag & m)
		{
			G=F ;
			// _MD = MatrixDomain<_Field>(F);
			VD = VectorDomain<_Field>(F);
			return;
		}


		void changeField(const _Field &F)
		{
			changeFieldSpecialised(const_cast<_Field&>(_field),
					       // const_cast<MatrixDomain<_Field>&>(_MD),
					       const_cast<VectorDomain<_Field>&>(_VD),
					       F,
					       typename FieldTraits<_Field>::categoryTag());
		}



	}; // end of class BlasMatrix


} // end of namespace LinBox

namespace LinBox
{ /* Blas Submatrix */
	/*! Dense Submatrix representation.
	 * @ingroup matrix
	 * A @ref BlasSubmatrix is a matrix of \p _Field::Element, with the structure of BLAS matrices.
	 * It is basically a read/write view on a vector of \p _Field::Element.
	 * In the Mother model, a @ref BlasSubmatrix is not allocated.
	 * <p>
	 * This matrix type conforms to the same interface as @ref BlasMatrix,
	 * except that you cannot resize it. It represents a submatrix of a dense
	 * matrix. Upon construction, one can freely manipulate the entries in the
	 * DenseSubmatrix, and the corresponding entries in the underlying
	 * @ref BlasMatrix will be modified.


	 */
	template <class _Matrix>
	class BlasSubmatrix {
	public :
		typedef typename _Matrix::Field           Field;
		typedef typename Field::Element         Element;    //!< Element type
		typedef typename _Matrix::Rep               Rep;    //!< Actually a <code>std::vector<Element></code> (or alike.)
		typedef BlasSubmatrix<typename _Matrix::Self_t>              Self_t;         //!< Self type
		typedef const BlasSubmatrix<typename _Matrix::constSelf_t>   constSelf_t;    //!< Self type (const)

		typedef typename Rep::pointer           pointer;    //!< pointer type to elements
		typedef const pointer             const_pointer;    //!< const pointer type
                typedef Self_t                    subMatrixType;    //!< Submatrix type
                typedef constSelf_t          constSubMatrixType;    //!< Submatrix type (const)
                typedef BlasMatrix<Field,Rep>             matrixType;    //!< matrix type
                typedef const BlasMatrix<Field,Rep>  constMatrixType;    //!< matrix type (const)
                typedef matrixType                     blasType;    //!< blas matrix type
                typedef BlasVector<Field,Rep>        vectorType;    //!< blas matrix type


	protected:
		_Matrix &_Mat;       //!< Parent BlasMatrix (ie encapsulated raw std::vector)
		size_t _row;                   //!< row dimension of Submatrix
		size_t _col;                   //!< col dimension of Submatrix
		size_t _r0;                    //!< upper left corner row of Submatrix in \p _Mat
		size_t _c0;                    //!< upper left corner row of Submatrix in \p _Mat
		size_t _stride ;               //!< number of columns in \p _Mat (or stride of \p _Mat)
		size_t _off;                   //!< offset in \p _Mat, precomputed \c (_row*_stride+_col)

		// applyDomain<matrixType>    _AD;
		applyDomain<constMatrixType>    _AD;
	public:
		VectorDomain<Field>    _VD; //!@bug NOT HERE

		//////////////////
		// CONSTRUCTORS //
		//////////////////


		/*  constructors */

		// /** NULL constructor.  */
		// BlasSubmatrix () ;

		/** Constructor from an existing @ref BlasMatrix and dimensions.
		 * \param M Pointer to @ref BlasMatrix of which to construct submatrix
		 * \param rowbeg Starting row
		 * \param colbeg Starting column
		 * \param Rowdim Row dimension
		 * \param Coldim Column dimension
		 */
		BlasSubmatrix (constMatrixType &M,
			       size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		BlasSubmatrix (matrixType &M,
			       size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);



		/** Constructor from an existing @ref BlasMatrix
		 * \param M Pointer to @ref BlasMatrix of which to construct submatrix
		 */
		BlasSubmatrix (constMatrixType &M);
		BlasSubmatrix (matrixType &M);

		//! @todo  BlasSub from (sub)Vector
		// BlasSubmatrix (const vectorType &V);


		/** Constructor from an existing submatrix and dimensions
		 * @param SM Constant reference to BlasSubmatrix from which to
		 *           construct submatrix
		 * @param rowbeg Starting row
		 * @param colbeg Starting column
		 * @param Rowdim Row dimension
		 * @param Coldim Column dimension
		 */
		BlasSubmatrix (constSelf_t  &SM,
				size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		BlasSubmatrix (Self_t  &SM,
				size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		/** Copy constructor.
		 * @param SM Submatrix to copy
		 */
		BlasSubmatrix (constSelf_t &SM);
		BlasSubmatrix (Self_t &SM);


		/*  Members  */

		/** Assignment operator.
		 * Assign the given submatrix to this one
		 * This is <i>only</i> renaming !
		 * There is no copy because BlasSubmatrix owns nothing.
		 * @param SM Submatrix to assign
		 * @return Reference to this submatrix
		 */
		BlasSubmatrix &operator = (const BlasSubmatrix<_Matrix> &SM);

		// function for repurposing Submatrices.
		BlasSubmatrix &submatrix(constSelf_t &SM,
				size_t rowbeg,
				size_t colbeg,
				size_t Rowdim,
				size_t Coldim);

		/// This is deep copy of the data, operator= is a shallow copy.
		template<class Matrix>
		BlasSubmatrix &copy( const Matrix & B);

		/// Swap contents.  Shapes must be the same.
		BlasSubmatrix &swap( Self_t & B);

		/// Overwrite with zeroes.
		BlasSubmatrix &zero();

		/// Overwrite with random elements.
		void random();

		template<class T>
		void random(const T&)
		{
			return random() ;
		}

		template<typename _Tp1, class _Rep2 = Rep>
		struct rebind ;

		//////////////////
		//  DIMENSIONS  //
		//////////////////

		/** Get the number of rows in the matrix
		 * @return Number of rows in matrix
		 */
		size_t rowdim () const;

		/** Get the number of columns in the matrix
		 * @return Number of columns in matrix
		 */
		size_t coldim () const ;

		/*! Get the stride of the matrix.
		 * @return stride of submatrix (number of cols of dense base matrix)
		 */
		size_t getStride() const;
		size_t stride() const { return getStride() ;}
		size_t offset() const { return _off ; }


		///////////////////
		//      I/O      //
		///////////////////

		/** Read the matrix from an input stream.
		 * @param file Input stream from which to read
		 * @bug reading a submatrix should not be allowed !!
		 */
		// template<class Field>
		std::istream& read (std::istream &file); // autodetect ?


		/** Write the matrix to an output stream.
		 * @param os Output stream to which to write
		 * @param f write in some format (@ref Tag::FileFormat::Format). Default is MM's.
		 */
		std::ostream &write (std::ostream &os,
				     LINBOX_enum (Tag::FileFormat) f = Tag::FileFormat::MatrixMarket )const;

		//////////////////
		//   ELEMENTS   //
		//////////////////

		/*! @internal
		 * Get read-only pointer to the matrix data.
		 */
		pointer getPointer() const ;

		const_pointer &getConstPointer() const ;


		/*! @internal
		 * Get write pointer to the matrix data.
		 * Data may be changed this way.
		 */
		pointer/* & */ getWritePointer() ;


		/** Set the entry at (i, j).
		 * @param i Row number, 0...rowdim () - 1
		 * @param j Column number 0...coldim () - 1
		 * @param a_ij Element to set
		 */
		void setEntry (size_t i, size_t j, const Element &a_ij) ;

		/** Get a writeable reference to an entry in the matrix.
		 * @param i Row index of entry
		 * @param j Column index of entry
		 * @return Reference to matrix entry
		 */
		Element &refEntry (size_t i, size_t j) ;

		/** Get a read-only individual entry from the matrix.
		 * @param i Row index
		 * @param j Column index
		 * @return Const reference to matrix entry
		 */
		const Element &getEntry (size_t i, size_t j) const ;

		/** Get an entry and store it in the given value.
		 * This form is more in the Linbox style and is provided for interface
		 * compatibility with other parts of the library
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @return Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const ;


		///////////////////
		//   ITERATORS   //
		///////////////////

		//! @name Forward declaration of Raw Iterators.
		//@{
		class Iterator  ;
		class ConstIterator ;

		class IndexedIterator ;
		class ConstIndexedIterator ;
		//@}


		/** @name typedef'd Row Iterators.
		 *\brief
		 * The row iterator gives the rows of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a row vector in dense format
		 * @{
		 */
		typedef typename matrixType::RowIterator            RowIterator;
		typedef typename matrixType::ConstRowIterator       ConstRowIterator;
		typedef typename matrixType::Row                    Row;
		typedef typename matrixType::ConstRow               ConstRow;
		//@} Row Iterators

		/** @name typedef'd Column Iterators.
		 *\brief
		 * The columns iterator gives the columns of the
		 * matrix in ascending order. Dereferencing the iterator yields
		 * a column vector in dense format
		 * @{
		 */
		typedef typename matrixType::ColIterator            ColIterator;
		typedef typename matrixType::ConstColIterator       ConstColIterator;
		typedef typename matrixType::Col                    Col;
		typedef typename matrixType::Column                 Column;
		typedef typename matrixType::ConstCol               ConstCol;
		//@} // Column Iterators



		RowIterator      rowBegin ();        //!< iterator to the begining of a row
		RowIterator      rowEnd ();          //!< iterator to the end of a row
		ConstRowIterator rowBegin () const;  //!< const iterator to the begining of a row
		ConstRowIterator rowEnd ()   const;  //!< const iterator to the end of a row

		ColIterator      colBegin ();
		ColIterator      colEnd ();
		ConstColIterator colBegin () const;
		ConstColIterator colEnd ()   const;

		Iterator      Begin ();
		Iterator      End ();
		ConstIterator Begin () const;
		ConstIterator End ()   const;



		IndexedIterator      IndexedBegin();
		IndexedIterator      IndexedEnd();
		ConstIndexedIterator IndexedBegin() const;
		ConstIndexedIterator IndexedEnd()   const;

		/*!  operator[].
		 * Retrieve a reference to a row
		 * @param i Row index
		 */
		Row      operator[] (size_t i) ;
		ConstRow operator[] (size_t i) const ;

		///////////////////
		//   BLACK BOX   //
		///////////////////


		//!@bug every vector we use here should have a stride/be blas vectors so it's not really templated by Vector1 Vector2 in general
		template <class Vector1, class Vector2>
		Vector1&  apply (Vector1& y, const Vector2& x) const
		{
			// std::cout << "prepare apply  subMatrix" << std::endl;
			// constSelf_t A(*this);
			// std::cout << "........................" << std::endl;
			// _AD.apply(Tag::Transpose::NoTrans,y,field().one,A,field().zero,x);
			_AD.apply(Tag::Transpose::NoTrans,y,field().one,*this,field().zero,x);
			// std::cout << "........done............" << std::endl;
			return y;
		}

		//! @bug use Matrix domain
		template <class Vector1, class Vector2>
		Vector1&  applyTranspose (Vector1& y, const Vector2& x) const
		{
			// std::cout << "prepare applyT subMatrix" << std::endl;
			// constSelf_t A(*this);
			// std::cout << "........................" << std::endl;
			// _AD.apply(Tag::Transpose::Trans,y,field().one,A,field().zero,x);
			_AD.apply(Tag::Transpose::Trans,y,field().one,*this,field().zero,x);
			// std::cout << "........done............" << std::endl;

			return y;
		}

		const Field& field() const { return _Mat.field() ;}
		// Field & field() { return _Mat.field(); }
	};

}

namespace LinBox
{ /* Triangular Matrix */
	//! Triangular BLAS matrix.
	template <class _Field, class _Storage >
	class TriangularBlasMatrix: public BlasMatrix<_Field,_Storage> {

	protected:

		LINBOX_enum (Tag::Shape)          _uplo; //!< upper or lower triangular
		LINBOX_enum (Tag::Diag)           _diag; //!< unit or non unit diagonal

	public:
		typedef _Field                       Field;
		typedef _Storage                         Rep;
		typedef typename Field::Element      Element;      //!< Element type
		typedef BlasMatrix<Field,Rep>           Father_t;
		typedef TriangularBlasMatrix<Field,Rep> Self_t;


		/*! Constructor for a new \c TriangularBlasMatrix.
		 * @param F
		 * @param m rows
		 * @param n cols
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const Field & F,
				      const size_t m, const size_t n,
				      LINBOX_enum (Tag::Shape) x=Tag::Shape::Upper,
				      LINBOX_enum (Tag::Diag) y= Tag::Diag::NonUnit) ;

		/*! Constructor from a \c BlasMatrix (copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (const BlasMatrix<Field,Rep>& A,
				      LINBOX_enum (Tag::Shape) x=Tag::Shape::Upper,
				      LINBOX_enum (Tag::Diag) y= Tag::Diag::NonUnit) ;

		/*! Constructor from a \c BlasMatrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		TriangularBlasMatrix (BlasMatrix<Field,Rep>& A,
				      LINBOX_enum (Tag::Shape) x=Tag::Shape::Upper,
				      LINBOX_enum (Tag::Diag) y= Tag::Diag::NonUnit) ;

		/*! Constructor from a \c TriangularBlasMatrix (copy).
		 * @param A matrix
		 */
		TriangularBlasMatrix (const TriangularBlasMatrix<Field,Rep>& A) ;

		/*! Generic constructor from a \c Matrix (no copy).
		 * @param A matrix
		 * @param y (non)unit diagonal
		 * @param x (upp/low)er matrix
		 */
		template<class Matrix>
		TriangularBlasMatrix (const Matrix& A,
				      LINBOX_enum (Tag::Shape) x=Tag::Shape::Upper,
				      LINBOX_enum (Tag::Diag) y= Tag::Diag::NonUnit) ;

		/// get the shape of the matrix (upper or lower)
		LINBOX_enum (Tag::Shape) getUpLo() const ;

		/// Is the diagonal implicitly unit ?
		LINBOX_enum (Tag::Diag) getDiag() const ;

	}; // end of class TriangularBlasMatrix

} // LinBox

#include "linbox/matrix/random-matrix.h"

namespace LinBox
{
	//! @bug does not work for submatrices.
	//! @todo b should be the random generator
	template<>
	template<>
	void BlasMatrix<PID_integer, Vector<PID_integer>::Dense >::random<unsigned>(const unsigned & b)
	{
		// std::cout << "randomized " <<  b << std::endl;
		RandomIntegerIter<false> R((unsigned)b);
		typedef RandomIntegerIter<false> IntRandIter ;
		typedef RandomDenseMatrix<IntRandIter, PID_integer > IntRand_t;
		IntRand_t Randomize(field(),R);
		Randomize.random(*this);

	}

} // LinBox

#include "blas-matrix.inl"
#include "blas-submatrix.inl"
#include "blas-triangularmatrix.inl"

#endif // __LINBOX_densematrix_blas_matrix_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
