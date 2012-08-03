/* linbox/matrix/matrix-domain.h
 * Copyright (C) 2002 Zhendong Wan, Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------------------------------
 * 2002-11-26  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Added detailed documentation, cleaned up the interface slightly, and added
 * support for matrix traits. Added read, write, neg, negin, axpy, and
 * matrix-vector and matrix-black box operations.
 * ------------------------------------------------------------
 *
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

#ifndef __LINBOX_matrix_domain_H
#define __LINBOX_matrix_domain_H

#include <iostream>
#include <vector>

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-domain.h"
//#include "linbox/matrix/blas-matrix.h"

namespace LinBox
{
	template<class Field> class BlasMatrix;
	template<class Field> class BlasSubmatrix;

	/** \brief For specializing matrix arithmetic
	 *
	 * This class defines matrix categories that allow us to specialize the matrix
	 * arithmetic in \ref MatrixDomain for different matrix representations. For
	 * example, a sparse matrix may have an efficient iterator over row vectors but
	 * not over column vectors. Therefore, an algorithm that tries to iterate over
	 * column vectors will run very slowly. Hence a specialization that avoids using
	 * column vectors is used instead.
	 */

	struct MatrixCategories {
		struct BlackboxTag { };
		struct RowMatrixTag : public virtual BlackboxTag { };
		struct ColMatrixTag : public virtual BlackboxTag { };
		struct RowColMatrixTag : public RowMatrixTag, public ColMatrixTag { };
	};

	template <class Matrix>
	struct MatrixTraits {
		typedef Matrix MatrixType;
		typedef typename MatrixCategories::BlackboxTag MatrixCategory;
	};

	/** \brief Helper class to allow specializations of certain matrix-vector products
	 *
	 * This class implements a method mulColSPD that multiplies a
	 * column-represented matrix by a dense vector
	 */
	template <class Field>
	class MVProductDomain {
	public:
		typedef typename Field::Element Element;

		MVProductDomain () {}

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense (const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const;
	};

	/** Class of matrix arithmetic functions.
	 *
	 * This class encapuslated matrix-matrix and matrix-vector operations, roughly
	 * equivalent to BLAS levels 2 and 3. The arithmetic methods are parameterized
	 * by matrix type so that they may be used the same way with sparse matrices,
	 * dense matrices, and dense submatrices. Except where otherwise noted, they
	 * require the matrix inputs to meet the \ref BlasMatrix archetype.
	 *
	 * These methods are specialized so that they can run efficiently with different
	 * matrix representations. If a matrix has an efficient row iterator, but not an
	 * efficient column iterator, a specialization that makes use of the former will
	 * be selected. This allows a great deal of flexibility when dealing with sparse
	 * matrix arithmetic.
	 *
	 * For all of the arithmetic operations that output matrices, it is assumed that
	 * the output matrix has an efficient row iterator. In typical use, the output
	 * matrix will be a \ref BlasMatrix or a \ref BlasSubmatrix, which has
	 * efficient row and column iterators. In particular, one should not perform
	 * these arithmetic operations outputting to a \ref SparseMatrixBase.
	 *
	 * There are other restrictions. See the method-specific documentation for more
	 * details.
	 */
	template <class Field>
	class MatrixDomain : public MVProductDomain<Field> {
	public:
		typedef size_t Index;
		typedef typename Field::Element Element;
		typedef Element Scalar;
		typedef std::vector<Element> Vector;
		// subvector
		typedef BlasMatrix<Field> Matrix;
		typedef BlasSubmatrix<Field> Submatrix;

		MatrixDomain () {/*std::cerr << "MD def cstor" << std::endl;*/ } 

		void init(const Field & F) { _field = &F; _VD.init(F); }

		/// Constructor.
		//! @param F field for MatrixDomain operations.
		MatrixDomain (const Field &F) :
			_field (&F), _VD (F)
		{ /*std::cerr << "MD cstor " << this << std::endl;*/ }

		/// Copy operator.
		MatrixDomain& operator= (const MatrixDomain& MD)
		{
			_field = MD._field;
			_VD = MD._VD;
			return *this;
		}

		/** Retrieve the underlying field.
		 * Return a reference to the field that this matrix domain
		 * object uses
		 * @returns reference to field
		 */
		//@{
		const Field &field () const
		{
			return *_field;
		}

		//@}

		/** Print matrix.
		 * @param  os  Output stream to which matrix is written.
		 * @param  A   Matrix.
		 * @returns reference to os.
		 */
		template <class Matrix_>
		inline std::ostream &write (std::ostream &os, const Matrix_ &A) const
		{
			return A.write (os);
		}

		/** Read matrix.
		 * @param  is  Input stream from which matrix is read.
		 * @param  A   Matrix.
		 * @returns reference to is.
		 */
		template <class Matrix_>
		inline std::istream &read (std::istream &is, Matrix_ &A) const
		{
			return A.read (is, _field);
		}

		/** Matrix copy
		 * B <- A.
		 * Copy the contents of the matrix B to the matrix A
		 *
		 * Both matrices must support the same iterators, row or column.
		 *
		 * @param B Matrix B
		 * @param A Matrix A
		 * @returns Reference to B
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1 &copy (Matrix1 &B, const Matrix2 &A) const
		{
			return copySpecialized (B, A,
						typename MatrixTraits<Matrix1>::MatrixCategory (),
						typename MatrixTraits<Matrix2>::MatrixCategory ());
		}

		/** Matrix equality.
		 * Test whether the matrices A and B are equal
		 * @param A Input vector
		 * @param B Input vector
		 * @returns true if and only if the matrices A and B are equal
		 */
		template <class Matrix1, class Matrix2>
		bool areEqual (const Matrix1 &A, const Matrix2 &B) const
		{
			return areEqualSpecialized (B, A,
						    typename MatrixTraits<Matrix1>::MatrixCategory (),
						    typename MatrixTraits<Matrix2>::MatrixCategory ());
		}

		/** Matrix equality with zero.
		 * @param A Input matrix
		 * @returns true if and only if the matrix A is zero
		 */
		template <class Matrix_>
		inline bool isZero (const Matrix_ &A) const
		{
			return isZeroSpecialized (A, typename MatrixTraits<Matrix_>::MatrixCategory ());
		}

		/** Matrix-matrix addition
		 * C <- A + B.
		 *
		 * Each of A, B, and C must support the same iterator, either row or
		 * column
		 *
		 * @param C Output matrix C
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to C
		 */
		template <class Matrix1, class Matrix2, class Matrix3>
		inline Matrix1& add (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
		{
			return addSpecialized (C, A, B,
					       typename MatrixTraits<Matrix1>::MatrixCategory (),
					       typename MatrixTraits<Matrix2>::MatrixCategory (),
					       typename MatrixTraits<Matrix3>::MatrixCategory ());
		}

		/** Matrix-matrix in-place addition
		 * A <- A + B.
		 *
		 * Each of A and B must support the same iterator, either row or column
		 *
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to A
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1& addin (Matrix1 &A, const Matrix2 &B) const
		{
			return addinSpecialized (A, B,
						 typename MatrixTraits<Matrix1>::MatrixCategory (),
						 typename MatrixTraits<Matrix2>::MatrixCategory ());
		}

		/** Matrix-matrix subtraction
		 * C <- A - B.
		 *
		 * Each of A, B, and C must support the same iterator, either row or
		 * column
		 *
		 * @param C Output matrix C
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to C
		 */
		template <class Matrix1, class Matrix2, class Matrix3>
		inline Matrix1 &sub (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
		{
			return subSpecialized (C, A, B,
					       typename MatrixTraits<Matrix1>::MatrixCategory (),
					       typename MatrixTraits<Matrix2>::MatrixCategory (),
					       typename MatrixTraits<Matrix3>::MatrixCategory ());
		}

		/** Matrix-matrix in-place subtraction
		 * A <- A - B.
		 *
		 * Each of A and B must support the same iterator, either row or column
		 *
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to A
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1 &subin (Matrix1 &A, const Matrix2 &B) const
		{
			return subinSpecialized (A, B,
						 typename MatrixTraits<Matrix1>::MatrixCategory (),
						 typename MatrixTraits<Matrix2>::MatrixCategory ());
		}

		/** Matrix negate
		 * B <- -A.
		 *
		 * Each of A and B must support the same iterator, either row or column
		 *
		 * @param B Output matrix B
		 * @param A Input matrix A
		 * @returns reference to B
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1 &neg (Matrix1 &B, const Matrix2 &A) const
		{
			return negSpecialized (B, A,
					       typename MatrixTraits<Matrix1>::MatrixCategory (),
					       typename MatrixTraits<Matrix2>::MatrixCategory ());
		}

		/** Matrix in-place negate
		 * A <- -A.
		 * @param A Input matrix A; result is stored here
		 */
		template <class Matrix_>
		inline Matrix_ &negin (Matrix_ &A) const
		{
			return neginSpecialized (A, typename MatrixTraits<Matrix_>::MatrixCategory ());
		}

		/** Matrix-matrix multiply
		 * C <- A * B.
		 *
		 * C must support both row and column iterators, and the vector
		 * representations must be dense. Examples of supported matrices are
		 * \ref BlasMatrix and \ref BlasSubmatrix.
		 *
		 * Either A or B, or both, may have limited iterators. However, either A
		 * must support row iterators or B must support column iterators. If
		 * both A and B lack support for an iterator (either row or column),
		 * then C must support the same type of iterator as A and B.
		 *
		 * @param C Output matrix C
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to C
		 */
		template <class Matrix1, class Matrix2, class Matrix3>
		inline Matrix1 &mul (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const
		{
			return mulSpecialized (C, A, B,
					       typename MatrixTraits<Matrix1>::MatrixCategory (),
					       typename MatrixTraits<Matrix2>::MatrixCategory (),
					       typename MatrixTraits<Matrix3>::MatrixCategory ());
		}

		/** Matrix-matrix in-place multiply on the left
		 * B <- A * B.
		 *
		 * B should support both row and column iterators, and must be dense. A
		 * must support row iterators.
		 *
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to B
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix2 &leftMulin (const Matrix1 &A, Matrix2 &B) const;

		/** Matrix-matrix in-place multiply on the right
		 * A <- A * B.
		 *
		 * A should support both row and column iterators, and must be dense. B
		 * must support column iterators.
		 *
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to A
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1 &rightMulin (Matrix1 &A, const Matrix2 &B) const;

		/** Matrix-matrix in-place multiply
		 * A <- A * B.
		 *
		 * This is an alias for \ref rightMulin
		 *
		 * @param A Input matrix A
		 * @param B Input matrix B
		 * @returns Reference to A
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1 &mulin (Matrix1 &A, const Matrix2 &B) const
		{
			return rightMulin (A, B);
		}

		/** Matrix-scalar multiply
		 * C <- B * a.
		 *
		 * Multiply B by the scalar element a and store the result in C. B and C
		 * must support the same iterators.
		 *
		 * @param C Output matrix C
		 * @param B Input matrix B
		 * @param a Input scalar a
		 * @returns Reference to C
		 */
		template <class Matrix1, class Matrix2>
		inline Matrix1 &mul (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const
		{
			return mulSpecialized (C, B, a,
					       typename MatrixTraits<Matrix1>::MatrixCategory (),
					       typename MatrixTraits<Matrix2>::MatrixCategory ());
		}

		/** Matrix-scalar in-place multiply
		 * B <- B * a.
		 *
		 * Multiply B by the scalar element a in-place.
		 *
		 * @param B Input matrix B
		 * @param a Input scalar a
		 * @returns Reference to B
		 */
		template <class Matrix_>
		inline Matrix_ &mulin (Matrix_ &B, const typename Field::Element &a) const
		{
			return mulinSpecialized (B, a, typename MatrixTraits<Matrix_>::MatrixCategory ());
		}

		/** Matrix-matrix in-place axpy
		 * Y <- Y + A*X.
		 *
		 * This function combines \ref mul and \ref add, eliminating the need
		 * for an additional temporary in expressions of the form $Y = Y +
		 * AX$. Only one row of additional storage is required. Y may have
		 * either efficient row iterators or efficient column iterators, and the
		 * same restrictions on A and X apply as in \ref mul.
		 *
		 * Note that no out-of-place axpy is provided, since it gives no
		 * benefit. One may just as easily multiply into the result and call
		 * \ref addin.
		 *
		 * @param Y Input matrix Y; result is stored here
		 * @param A Input matrix A
		 * @param X Input matrix X
		 */
		template <class Matrix1, class Matrix2, class Matrix3>
		inline Matrix1 &axpyin (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const
		{
			return axpyinSpecialized (Y, A, X,
						  typename MatrixTraits<Matrix1>::MatrixCategory (),
						  typename MatrixTraits<Matrix2>::MatrixCategory (),
						  typename MatrixTraits<Matrix3>::MatrixCategory ());
		}

		//! Y <- AX-Y
		template <class Matrix1, class Matrix2, class Matrix3>
		inline Matrix1 &axmyin (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const
		{
			negin(Y);
			axpyin(Y,A,X);
			return Y;
		}


		/*!  General matrix multiply
		 * \f$ D \gets \alpha A B + \beta C\f$.
		 * @todo not efficient...
		 */
		template <class Matrix1, class Matrix2, class Matrix3>
		inline Matrix1 &muladd (Matrix1                       & D,
					const typename Field::Element & beta,
					const Matrix1                 & C,
					const typename Field::Element & alpha,
					const Matrix2                 & A,
					const Matrix3                 & B) const
		{
			mul(D,A,B); // D = AB
			mulin(D,alpha); // D = alpha D
			Matrix1 CC(C);
			mulin(CC,beta);   // C = beta C
			addin(D,CC);   // D = D+C
			return D;
		}

		/*! @todo Need documentation of these methods */
		//@{
		template<class Matrix1, class Matrix2>
		Matrix1 &pow_apply (Matrix1 &M1, const Matrix2 &M2, unsigned long int k) const;

		template<class Matrix1, class Matrix2>
		Matrix1 &pow_horn (Matrix1 &M1, const Matrix2 &M2, unsigned long int k) const;
		//@}


		/*! @name Matrix-vector arithmetic operations
		 * These operations take a matrix satisfying the \ref DenseMatrix
		 * archetype and LinBox vectors as inputs. They involve matrix-vector
		 * product and matrix-vector AXPY
		 */
		//@{
		/** Matrix-vector multiply
		 * w <- A * v.
		 *
		 * The vectors v and w must be of the same representation (dense, sparse
		 * sequence, sparse associative, or sparse parallel), but they may be of
		 * different types. The matrix A may have any representation.
		 *
		 * @param w Output vector w
		 * @param A Input matrix A
		 * @param v Input vector v
		 * @returns Reference to w
		 */
		template <class Vector1, class Matrix_, class Vector2>
		inline Vector1 &vectorMul (Vector1 &w, const Matrix_ &A, const Vector2 &v) const
		{
			return mulSpecialized (w, A, v, typename MatrixTraits<Matrix_>::MatrixCategory ());
		}

		/** Matrix-vector in-place axpy
		 * \f$y \gets y + A x\f$.
		 *
		 * This function eliminates the requirement for temporary storage when
		 * one is computing an expression of the form given above.
		 *
		 * The vectors y and x must be of the same representation (dense, sparse
		 * sequence, sparse associative, or sparse parallel), but they may be of
		 * different types. The matrix A may have any representation.
		 *
		 * Note that out-of-place axpy is not provided since it provides no
		 * benefit -- one can use mul and then addin to exactly the same effect,
		 * with no additional storage or performance cost.
		 *
		 * @param y Input vector y; result is stored here
		 * @param A Input matrix A
		 * @param x Input vector x
		 */
		template <class Vector1, class Matrix_, class Vector2>
		inline Vector1 &vectorAxpyin (Vector1 &y, const Matrix_ &A, const Vector2 &x) const
		{
			return axpyinSpecialized (y, A, x, typename MatrixTraits<Matrix_>::MatrixCategory ());
		}
		//@}

		/*! @name Matrix-black box arithmetic operations
		 * These operations mimic the matrix-matrix arithmetic operations above,
		 * but one of the parameters is a \ref BlackboxArchetype.
		 */
		//@{
		/** Matrix-black box left-multiply
		 * C <- A * B.
		 *
		 * Both C and B must support column iterators
		 *
		 * @param C Output matrix
		 * @param A Black box for A
		 * @param B Matrix B
		 */
		template <class Matrix1, class Blackbox, class Matrix2>
		inline Matrix1 &blackboxMulLeft (Matrix1 &C, const Blackbox &A, const Matrix2 &B) const;

		/** Matrix-black box right-multiply
		 * C <- A * B.
		 *
		 * Both C and A must support row iterators
		 *
		 * @param C Output matrix
		 * @param A Matrix A
		 * @param B Black box for B
		 */
		template <class Matrix1, class Matrix2, class Blackbox>
		inline Matrix1 &blackboxMulRight (Matrix1 &C, const Matrix2 &A, const Blackbox &B) const;
		//@}

		/*! @name Matrix permutations
		 * @brief
		 * These operations permute the rows or columns of a matrix based on
		 * the given permutation. They are intended for use with Gauss-Jordan
		 * elimination
		 */
		//@{
		/// Transposition.
		typedef std::pair<unsigned int, unsigned int> Transposition;
		/** Permutation.
		 *
		 * A permutation is represented as a vector of pairs, each
		 * pair representing a transposition.
		 */
		typedef std::vector<Transposition> Permutation;


		/** Permute the rows of the given matrix.
		 *
		 * @param A Output matrix
		 * @param P_start Start of permutation
		 * @param P_end End of permutation
		 * @returns Reference to A
		 */
		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteRows (Matrix_   &A,
					    Iterator  P_start,
					    Iterator  P_end) const
		{
			return permuteRowsSpecialized (A, P_start, P_end,
						       typename MatrixTraits<Matrix_>::MatrixCategory ());
		}

		/** Permute the columns of the given matrix.
		 *
		 * @param A Output matrix
		 * @param P_start Start of permutation
		 * @param P_end End of permutation
		 * @returns Reference to A
		 */
		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteColumns (Matrix_   &A,
					       Iterator  P_start,
					       Iterator  P_end) const
		{
			return permuteColsSpecialized (A, P_start, P_end,
						       typename MatrixTraits<Matrix_>::MatrixCategory ());
		}
		//@}

	protected:

		// Specialized function implementations
		template <class Matrix1, class Matrix2> Matrix1 &copyRow (Matrix1 &B, const Matrix2 &A) const;
		template <class Matrix1, class Matrix2> Matrix1 &copyCol (Matrix1 &B, const Matrix2 &A) const;

		template <class Matrix1, class Matrix2>
		inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
						 MatrixCategories::RowMatrixTag,
						 MatrixCategories::RowMatrixTag) const
		{
			return copyRow (B, A);
		}
		template <class Matrix1, class Matrix2>
		inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
						 MatrixCategories::ColMatrixTag,
						 MatrixCategories::ColMatrixTag) const
		{
			return copyCol (B, A);
		}
		template <class Matrix1, class Matrix2>
		inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
						 MatrixCategories::RowColMatrixTag,
						 MatrixCategories::RowColMatrixTag) const
		{
			return copyRow (B, A);
		}

		template <class Matrix1, class Matrix2> bool areEqualRow (const Matrix1 &A, const Matrix2 &B) const;
		template <class Matrix1, class Matrix2> bool areEqualCol (const Matrix1 &A, const Matrix2 &B) const;

		template <class Matrix1, class Matrix2>
		inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
						 MatrixCategories::RowMatrixTag,
						 MatrixCategories::RowMatrixTag) const
		{
			return areEqualRow (A, B);
		}
		template <class Matrix1, class Matrix2>
		inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
						 MatrixCategories::ColMatrixTag,
						 MatrixCategories::ColMatrixTag) const
		{
			return areEqualCol (A, B);
		}
		template <class Matrix1, class Matrix2>
		inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
						 MatrixCategories::RowColMatrixTag,
						 MatrixCategories::RowColMatrixTag) const
		{
			return areEqualRow (A, B);
		}

		template <class Matrix_> bool isZeroRow (const Matrix_ &v) const;
		template <class Matrix_> bool isZeroCol (const Matrix_ &v) const;

		template <class Matrix_>
		bool isZeroSpecialized (const Matrix_ &A, MatrixCategories::RowMatrixTag) const
		{
			return isZeroRow (A);
		}
		template <class Matrix_>
		bool isZeroSpecialized (const Matrix_ &A, MatrixCategories::ColMatrixTag) const
		{
			return isZeroCol (A);
		}
		template <class Matrix_>
		bool isZeroSpecialized (const Matrix_ &A, MatrixCategories::RowColMatrixTag) const
		{
			return isZeroRow (A);
		}

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& addRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& addCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& addSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{
			return addRow (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& addSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return addCol (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& addSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag) const
		{
			return addRow (C, A, B);
		}

		template <class Matrix1, class Matrix2> Matrix1& addinRow (Matrix1 &A, const Matrix2 &B) const;
		template <class Matrix1, class Matrix2> Matrix1& addinCol (Matrix1 &A, const Matrix2 &B) const;

		template <class Matrix1, class Matrix2>
		inline Matrix1& addinSpecialized (Matrix1 &A, const Matrix2 &B,
						  MatrixCategories::RowMatrixTag,
						  MatrixCategories::RowMatrixTag) const
		{
			return addinRow (A, B);
		}
		template <class Matrix1, class Matrix2>
		inline Matrix1& addinSpecialized (Matrix1 &A, const Matrix2 &B,
						  MatrixCategories::ColMatrixTag,
						  MatrixCategories::ColMatrixTag) const
		{
			return addinCol (A, B);
		}
		template <class Matrix1, class Matrix2>
		inline Matrix1& addinSpecialized (Matrix1 &A, const Matrix2 &B,
						  MatrixCategories::RowColMatrixTag,
						  MatrixCategories::RowColMatrixTag) const
		{
			return addinRow (A, B);
		}

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& subRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& subCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& subSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{
			return subRow (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& subSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return subCol (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1& subSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag) const
		{
			return subRow (C, A, B);
		}

		template <class Matrix1, class Matrix2> Matrix1& subinRow (Matrix1 &A, const Matrix2 &B) const;
		template <class Matrix1, class Matrix2> Matrix1& subinCol (Matrix1 &A, const Matrix2 &B) const;

		template <class Matrix1, class Matrix2>
		Matrix1& subinSpecialized (Matrix1 &A, const Matrix2 &B,
					   MatrixCategories::RowMatrixTag,
					   MatrixCategories::RowMatrixTag) const
		{
			return subinRow (A, B);
		}
		template <class Matrix1, class Matrix2>
		Matrix1& subinSpecialized (Matrix1 &A, const Matrix2 &B,
					   MatrixCategories::ColMatrixTag,
					   MatrixCategories::ColMatrixTag) const
		{
			return subinCol (A, B);
		}
		template <class Matrix1, class Matrix2>
		Matrix1& subinSpecialized (Matrix1 &A, const Matrix2 &B,
					   MatrixCategories::RowColMatrixTag,
					   MatrixCategories::RowColMatrixTag) const
		{
			return subinRow (A, B);
		}

		template <class Matrix1, class Matrix2> Matrix1& negRow (Matrix1 &A, const Matrix2 &B) const;
		template <class Matrix1, class Matrix2> Matrix1& negCol (Matrix1 &A, const Matrix2 &B) const;

		template <class Matrix1, class Matrix2>
		inline Matrix1& negSpecialized (Matrix1 &A, const Matrix2 &B,
						MatrixCategories::RowMatrixTag,
						MatrixCategories::RowMatrixTag) const
		{
			return negRow (A, B);
		}
		template <class Matrix1, class Matrix2>
		inline Matrix1& negSpecialized (Matrix1 &A, const Matrix2 &B,
						MatrixCategories::ColMatrixTag,
						MatrixCategories::ColMatrixTag) const
		{
			return negCol (A, B);
		}
		template <class Matrix1, class Matrix2>
		inline Matrix1& negSpecialized (Matrix1 &A, const Matrix2 &B,
						MatrixCategories::RowColMatrixTag,
						MatrixCategories::RowColMatrixTag) const
		{
			return negRow (A, B);
		}

		template <class Matrix_> Matrix_ &neginRow (Matrix_ &A) const;
		template <class Matrix_> Matrix_ &neginCol (Matrix_ &A) const;

		template <class Matrix_>
		Matrix_ &neginSpecialized (Matrix_ &A, MatrixCategories::RowMatrixTag) const
		{
			return neginRow (A);
		}
		template <class Matrix_>
		Matrix_ &neginSpecialized (Matrix_ &A, MatrixCategories::ColMatrixTag) const
		{
			return neginCol (A);
		}
		template <class Matrix_>
		Matrix_ &neginSpecialized (Matrix_ &A, MatrixCategories::RowColMatrixTag) const
		{
			return neginRow (A);
		}

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulRowRowCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulColRowCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulRowRowRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulColColCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::BlackboxTag) const
		{
			return blackboxMulRight(C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::BlackboxTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return blackboxMulLeft(C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulRowRowCol (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulColRowCol (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulRowRowCol (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{
			return mulRowRowRow (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulColColCol (C, A, B);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag) const
		{
			return mulRowRowCol (C, A, B);
		}

		template <class Matrix1, class Matrix2>
		Matrix1 &mulRow (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const;
		template <class Matrix1, class Matrix2>
		Matrix1 &mulCol (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const;

		template <class Matrix1, class Matrix2>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a,
					 MatrixCategories::RowMatrixTag,
					 MatrixCategories::RowMatrixTag) const
		{
			return mulRow (C, B, a);
		}
		template <class Matrix1, class Matrix2>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a,
					 MatrixCategories::ColMatrixTag,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulCol (C, B, a);
		}
		template <class Matrix1, class Matrix2>
		Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag) const
		{
			return mulRow (C, B, a);
		}

		template <class Matrix_> Matrix_ &mulinRow (Matrix_ &B, const typename Field::Element &a) const;
		template <class Matrix_> Matrix_ &mulinCol (Matrix_ &B, const typename Field::Element &a) const;

		template <class Matrix_>
		Matrix_ &mulinSpecialized (Matrix_ &B, const typename Field::Element &a,
					  MatrixCategories::RowMatrixTag) const
		{
			return mulinRow (B, a);
		}
		template <class Matrix_>
		Matrix_ &mulinSpecialized (Matrix_ &B, const typename Field::Element &a,
					  MatrixCategories::ColMatrixTag) const
		{
			return mulinCol (B, a);
		}
		template <class Matrix_>
		Matrix_ &mulinSpecialized (Matrix_ &B, const typename Field::Element &a,
					  MatrixCategories::RowColMatrixTag) const
		{
			return mulinRow (B, a);
		}

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinRowRowCol (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinColRowCol (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinRowRowRow (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const;
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinColColCol (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X) const;

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
					    MatrixCategories::RowMatrixTag,
					    MatrixCategories::RowMatrixTag,
					    MatrixCategories::ColMatrixTag) const
		{
			return axpyinRowRowCol (Y, A, X);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
					    MatrixCategories::ColMatrixTag,
					    MatrixCategories::RowMatrixTag,
					    MatrixCategories::ColMatrixTag) const
		{
			return axpyinColRowCol (Y, A, X);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
					    MatrixCategories::RowColMatrixTag,
					    MatrixCategories::RowMatrixTag,
					    MatrixCategories::ColMatrixTag) const
		{
			return axpyinRowRowCol (Y, A, X);
		}
		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
					    MatrixCategories::RowMatrixTag,
					    MatrixCategories::RowMatrixTag,
					    MatrixCategories::RowMatrixTag) const
		{
			return axpyinRowRowRow (Y, A, X);
		}

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
					    MatrixCategories::ColMatrixTag,
					    MatrixCategories::ColMatrixTag,
					    MatrixCategories::ColMatrixTag) const
		{
			return axpyinColColCol (Y, A, X);
		}

		template <class Matrix1, class Matrix2, class Matrix3>
		Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
					    MatrixCategories::RowColMatrixTag,
					    MatrixCategories::RowColMatrixTag,
					    MatrixCategories::RowColMatrixTag) const
		{
			return axpyinRowRowCol (Y, A, X);
		}

		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulRowSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulRowSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulRowSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulRowSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::SparseParallelVectorTag) const;

		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulColSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::DenseVectorTag,
					    VectorCategories::DenseVectorTag) const
		{
			return this->mulColDense (_VD, w, A, v);
		}
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulColSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::DenseVectorTag,
					    VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulColSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::DenseVectorTag,
					    VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulColSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					    VectorCategories::DenseVectorTag,
					    VectorCategories::SparseParallelVectorTag) const;

		template <class Vector1, class Matrix_, class Vector2>
		inline Vector1 &mulColSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
						   VectorCategories::GenericVectorTag,
						   VectorCategories::GenericVectorTag) const
		{
			typename LinBox::Vector<Field>::Dense y;

			VectorWrapper::ensureDim (y, w.size ());

			VectorWrapper::ensureDim (y, w.size ());

			vectorMul (y, A, v);
			_VD.copy (w, y);

			return w;
		}

		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					 MatrixCategories::RowMatrixTag) const
		{
			return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ());
		}
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					 MatrixCategories::ColMatrixTag) const
		{
			return mulColSpecialized (w, A, v,
						  typename VectorTraits<Vector1>::VectorCategory (),
						  typename VectorTraits<Vector2>::VectorCategory ());
		}
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &mulSpecialized (Vector1 &w, const Matrix_ &A, const Vector2 &v,
					 MatrixCategories::RowColMatrixTag) const
		{
			return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ());
		}

		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::SparseParallelVectorTag) const;

		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					       VectorCategories::SparseParallelVectorTag) const;

		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					    MatrixCategories::RowMatrixTag) const
		{
			return axpyinRowSpecialized (y, A, x, typename VectorTraits<Vector1>::VectorCategory ());
		}
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					    MatrixCategories::ColMatrixTag) const
		{
			return axpyinColSpecialized (y, A, x, typename VectorTraits<Vector1>::VectorCategory ());
		}
		template <class Vector1, class Matrix_, class Vector2>
		Vector1 &axpyinSpecialized (Vector1 &y, const Matrix_ &A, const Vector2 &x,
					    MatrixCategories::RowColMatrixTag) const
		{
			return axpyinRowSpecialized (y, A, x, typename VectorTraits<Vector1>::VectorCategory ());
		}

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteRowsByRow (Matrix_   &A,
						 Iterator  P_start,
						 Iterator  P_end) const;

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteRowsByCol (Matrix_   &A,
						 Iterator  P_start,
						 Iterator  P_end) const;

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteRowsSpecialized (Matrix_   &A,
						       Iterator  P_start,
						       Iterator  P_end,
						       MatrixCategories::RowColMatrixTag) const
		{
			return permuteRowsByCol (A, P_start, P_end);
		}

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteRowsSpecialized (Matrix_   &A,
						       Iterator  P_start,
						       Iterator  P_end,
						       MatrixCategories::RowMatrixTag) const
		{
			return permuteRowsByRow (A, P_start, P_end);
		}

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteRowsSpecialized (Matrix_   &A,
						       Iterator  P_start,
						       Iterator  P_end,
						       MatrixCategories::ColMatrixTag) const
		{
			return permuteRowsByCol (A, P_start, P_end);
		}

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteColsByRow (Matrix_   &A,
						 Iterator  P_start,
						 Iterator  P_end) const;

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteColsByCol (Matrix_   &A,
						 Iterator  P_start,
						 Iterator  P_end) const;

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteColsSpecialized (Matrix_   &A,
						       Iterator  P_start,
						       Iterator  P_end,
						       MatrixCategories::RowColMatrixTag) const
		{
			return permuteColsByRow (A, P_start, P_end);
		}

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteColsSpecialized (Matrix_   &A,
						       Iterator  P_start,
						       Iterator  P_end,
						       MatrixCategories::RowMatrixTag) const
		{
			return permuteColsByRow (A, P_start, P_end);
		}

		template <class Matrix_, class Iterator>
		inline Matrix_ &permuteColsSpecialized (Matrix_   &A,
						       Iterator  P_start,
						       Iterator  P_end,
						       MatrixCategories::ColMatrixTag) const
		{
			return permuteColsByCol (A, P_start, P_end);
		}

		const Field         *_field;
		VectorDomain<Field>  _VD;
	}; //MatrixDomain

}

#include "linbox/matrix/matrix-domain.inl"

#endif // __LINBOX_matrix_domain_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

