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
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_domain_H
#define __LINBOX_matrix_domain_H

#include <iostream>

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-domain.h"

namespace LinBox
{

/** \brief For specializing matrix arithmetic
 *
 * This class defines matrix categories that allow us to specialize the matrix
 * arithmetic in \ref MatrixDomain for different matrix representations. For
 * example, a sparse matrix may have an efficient iterator over row vectors but
 * not over column vectors. Therefore, an algorithm that tries to iterate over
 * column vectors will run very slowly. Hence a specialization that avoids using
 * column vectors is used instead.
 */

struct MatrixCategories 
{
	struct BlackboxTag { };
	struct RowMatrixTag : public virtual BlackboxTag { };
	struct ColMatrixTag : public virtual BlackboxTag { };
	struct RowColMatrixTag : public RowMatrixTag, public ColMatrixTag { };
};

template <class Matrix> struct MatrixTraits
{
	typedef Matrix MatrixType;
	typedef typename Matrix::MatrixCategory MatrixCategory;
};

/** \brief Helper class to allow specializations of certain matrix-vector products
 *
 * This class implements a method mulColSPD that multiplies a
 * column-represented matrix by a dense vector
 */

template <class Field>
class MVProductDomain
{
    public:
	typedef typename Field::Element Element;

	MVProductDomain () {}

    protected:
	template <class Vector1, class Matrix, class Vector2>
	inline Vector1 &mulColDense (const VectorDomain<Field> &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const;
};

/** \brief Class of matrix arithmetic functions
 *
 * This class encapuslated matrix-matrix and matrix-vector operations, roughly
 * equivalent to BLAS levels 2 and 3. The arithmetic methods are parameterized
 * by matrix type so that they may be used the same way with sparse matrices,
 * dense matrices, and dense submatrices. Except where otherwise noted, they
 * require the matrix inputs to meet the \ref DenseMatrix archetype.
 *
 * These methods are specialized so that they can run efficiently with different
 * matrix representations. If a matrix has an efficient row iterator, but not an
 * efficient column iterator, a specialization that makes use of the former will
 * be selected. This allows a great deal of flexibility when dealing with sparse
 * matrix arithmetic.
 *
 * For all of the arithmetic operations that output matrices, it is assumed that
 * the output matrix has an efficient row iterator. In typical use, the output
 * matrix will be a \ref DenseMatrixBase or a \ref DenseSubmatrix, which has
 * efficient row and column iterators. In particular, one should not perform
 * these arithmetic operations outputting to a \ref SparseMatrixBase.
 *
 * There are other restrictions. See the method-specific documentation for more
 * details.
 */

template <class Field>
class MatrixDomain : public MVProductDomain<Field>
{
    public:

	///
	MatrixDomain (const Field &F) : _F (F), _VD (F) {}

	MatrixDomain& operator= (const MatrixDomain& MD)
	{ _F = MD._F; 
	  _VD = MD._VD; 
	  return *this; 
	}

	/** Retrieve the underlying field
	 * Return a reference to the field that this matrix domain
	 * object uses
	 * @returns reference to field
	 */
	const Field &field () const
		{ return _F; }
	Field &field () 
		{ return _F; }

	/** Print matrix.
	 * @param  os  Output stream to which matrix is written.
	 * @param  A   Matrix.
	 * @returns reference to os.
	 */
	template <class Matrix>
	inline std::ostream &write (std::ostream &os, const Matrix &A) const
		{ return A.write (os, _F); }

	/** Read matrix
	 * @param  is  Input stream from which matrix is read.
	 * @param  A   Matrix.
	 * @returns reference to is.
	 */
	template <class Matrix>
	inline std::istream &read (std::istream &is, Matrix &A) const
		{ return A.read (is, _F); }

	/** Matrix copy
	 * B <- A
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
		{ return copySpecialized (B, A,
					  typename MatrixTraits<Matrix1>::MatrixCategory (),
					  typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix equality
	 * Test whether the matrices A and B are equal
	 * @param A Input vector
	 * @param B Input vector
	 * @returns true if and only if the matrices A and B are equal
	 */
	template <class Matrix1, class Matrix2>
	bool areEqual (const Matrix1 &A, const Matrix2 &B) const
		{ return areEqualSpecialized (B, A,
					      typename MatrixTraits<Matrix1>::MatrixCategory (),
					      typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix equality with zero
	 * @param A Input matrix
	 * @returns true if and only if the matrix A is zero
	 */
	template <class Matrix>
	inline bool isZero (const Matrix &A) const
		{ return isZeroSpecialized (A, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Matrix-matrix addition
	 * C <- A + B
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
		{ return addSpecialized (C, A, B,
					 typename MatrixTraits<Matrix1>::MatrixCategory (),
					 typename MatrixTraits<Matrix2>::MatrixCategory (),
					 typename MatrixTraits<Matrix3>::MatrixCategory ()); }
    
	/** Matrix-matrix in-place addition
	 * A <- A + B
	 *
	 * Each of A and B must support the same iterator, either row or column
	 *
	 * @param A Input matrix A
	 * @param B Input matrix B
	 * @returns Reference to A
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix1& addin (Matrix1 &A, const Matrix2 &B) const
		{ return addinSpecialized (A, B,
					   typename MatrixTraits<Matrix1>::MatrixCategory (),
					   typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix-matrix subtraction
	 * C <- A - B
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
		{ return subSpecialized (C, A, B,
					 typename MatrixTraits<Matrix1>::MatrixCategory (),
					 typename MatrixTraits<Matrix2>::MatrixCategory (),
					 typename MatrixTraits<Matrix3>::MatrixCategory ()); }

	/** Matrix-matrix in-place subtraction
	 * A <- A - B
	 *
	 * Each of A and B must support the same iterator, either row or column
	 *
	 * @param A Input matrix A
	 * @param B Input matrix B
	 * @returns Reference to A
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix1 &subin (Matrix1 &A, const Matrix2 &B) const
		{ return subinSpecialized (A, B,
					   typename MatrixTraits<Matrix1>::MatrixCategory (),
					   typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix negate
	 * B <- -A
	 *
	 * Each of A and B must support the same iterator, either row or column
	 *
	 * @param B Output matrix B
	 * @param A Input matrix A
	 * @returns reference to B
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix1 &neg (Matrix1 &B, const Matrix2 &A) const
		{ return negSpecialized (B, A,
					 typename MatrixTraits<Matrix1>::MatrixCategory (),
					 typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix in-place negate
	 * A <- -A
	 * @param A Input matrix A; result is stored here
	 */
	template <class Matrix>
	inline Matrix &negin (Matrix &A) const
		{ return neginSpecialized (A, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Matrix-matrix multiply
	 * C <- A * B
	 *
	 * C must support both row and column iterators, and the vector
	 * representations must be dense. Examples of supported matrices are
	 * \ref DenseMatrixBase and \ref DenseSubmatrix.
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
		{ return mulSpecialized (C, A, B,
					 typename MatrixTraits<Matrix1>::MatrixCategory (),
					 typename MatrixTraits<Matrix2>::MatrixCategory (),
					 typename MatrixTraits<Matrix3>::MatrixCategory ()); }

	/** Matrix-matrix in-place multiply on the left
	 * B <- A * B
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
	 * A <- A * B
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
	 * A <- A * B
	 *
	 * This is an alias for \ref rightMulin
	 *
	 * @param A Input matrix A
	 * @param B Input matrix B
	 * @returns Reference to A
	 */
	template <class Matrix1, class Matrix2>
	inline Matrix1 &mulin (Matrix1 &A, const Matrix2 &B) const
		{ return rightMulin (A, B); }

	/** Matrix-scalar multiply
	 * C <- B * a
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
		{ return mulSpecialized (C, B, a,
					 typename MatrixTraits<Matrix1>::MatrixCategory (),
					 typename MatrixTraits<Matrix2>::MatrixCategory ()); }

	/** Matrix-scalar in-place multiply
	 * B <- B * a
	 *
	 * Multiply B by the scalar element a in-place.
	 *
	 * @param B Input matrix B
	 * @param a Input scalar a
	 * @returns Reference to B
	 */
	template <class Matrix>
	inline Matrix &mulin (Matrix &B, const typename Field::Element &a) const
		{ return mulinSpecialized (B, a, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Matrix-matrix in-place axpy
	 * Y <- Y + A*X
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
		{ return axpyinSpecialized (Y, A, X,
					    typename MatrixTraits<Matrix1>::MatrixCategory (),
					    typename MatrixTraits<Matrix2>::MatrixCategory (),
					    typename MatrixTraits<Matrix3>::MatrixCategory ()); }

	/* FIXME: Need documentation of these methods */
	template<class Matrix1, class Matrix2>
	Matrix1 &pow_apply (Matrix1 &M1, const Matrix2 &M2, unsigned long int k) const;

	template<class Matrix1, class Matrix2>
	Matrix1 &pow_horn (Matrix1 &M1, const Matrix2 &M2, unsigned long int k) const;


	/* @name Matrix-vector arithmetic operations
	 * These operations take a matrix satisfying the \ref DenseMatrix
	 * archetype and LinBox vectors as inputs. They involve matrix-vector
	 * product and matrix-vector AXPY
	 */

	/** Matrix-vector multiply
	 * w <- A * v
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
	template <class Vector1, class Matrix, class Vector2>
	inline Vector1 &vectorMul (Vector1 &w, const Matrix &A, const Vector2 &v) const
		{ return mulSpecialized (w, A, v, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Matrix-vector in-place axpy
	 * y <- y + A*x
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
	template <class Vector1, class Matrix, class Vector2>
	inline Vector1 &vectorAxpyin (Vector1 &y, const Matrix &A, const Vector2 &x) const
		{ return axpyinSpecialized (y, A, x, typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/*? @name Matrix-black box arithmetic operations
	 * These operations mimic the matrix-matrix arithmetic operations above,
	 * but one of the parameters is a \ref BlackboxArchetype.
	 */

	/** Matrix-black box left-multiply
	 * C <- A * B
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
	 * C <- A * B
	 *
	 * Both C and A must support row iterators
	 *
	 * @param C Output matrix
	 * @param A Matrix A
	 * @param B Black box for B
	 */
	template <class Matrix1, class Matrix2, class Blackbox>
	inline Matrix1 &blackboxMulRight (Matrix1 &C, const Matrix2 &A, const Blackbox &B) const;

	/*? @name Matrix permutations
	 * These operations permute the rows or columns of a matrix based on
	 * the given permutation. They are intended for use with Gauss-Jordan
	 * elimination
	 */

	/** Permutation
	 *
	 * A permutation is represented as a vector of pairs, each
	 * pair representing a transposition.
	 */
	typedef std::pair<unsigned int, unsigned int> Transposition;
	typedef std::vector<Transposition> Permutation;



	/** Permute the rows of the given matrix
	 *
	 * @param A Output matrix
	 * @param P_start Start of permutation
	 * @param P_end End of permutation
	 * @returns Reference to A
	 */
	template <class Matrix, class Iterator>
	inline Matrix &permuteRows (Matrix   &A,
				    Iterator  P_start,
				    Iterator  P_end) const
		{ return permuteRowsSpecialized (A, P_start, P_end,
						 typename MatrixTraits<Matrix>::MatrixCategory ()); }

	/** Permute the columns of the given matrix
	 *
	 * @param A Output matrix
	 * @param P_start Start of permutation
	 * @param P_end End of permutation
	 * @returns Reference to A
	 */
	template <class Matrix, class Iterator>
	inline Matrix &permuteColumns (Matrix   &A,
				       Iterator  P_start,
				       Iterator  P_end) const
		{ return permuteColsSpecialized (A, P_start, P_end,
						 typename MatrixTraits<Matrix>::MatrixCategory ()); }

    private:

	// Specialized function implementations
	template <class Matrix1, class Matrix2> Matrix1 &copyRow (Matrix1 &B, const Matrix2 &A) const;
	template <class Matrix1, class Matrix2> Matrix1 &copyCol (Matrix1 &B, const Matrix2 &A) const;

	template <class Matrix1, class Matrix2>
	inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag) const
		{ return copyRow (B, A); }
	template <class Matrix1, class Matrix2>
	inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return copyCol (B, A); }
	template <class Matrix1, class Matrix2>
	inline Matrix1 &copySpecialized (Matrix1 &B, const Matrix2 &A,
					 MatrixCategories::RowColMatrixTag,
					 MatrixCategories::RowColMatrixTag) const
		{ return copyRow (B, A); }

	template <class Matrix1, class Matrix2> bool areEqualRow (const Matrix1 &A, const Matrix2 &B) const;
	template <class Matrix1, class Matrix2> bool areEqualCol (const Matrix1 &A, const Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
				  MatrixCategories::RowMatrixTag,
				  MatrixCategories::RowMatrixTag) const
		{ return areEqualRow (A, B); }
	template <class Matrix1, class Matrix2>
	inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
				  MatrixCategories::ColMatrixTag,
				  MatrixCategories::ColMatrixTag) const
		{ return areEqualCol (A, B); }
	template <class Matrix1, class Matrix2>
	inline bool areEqualSpecialized (const Matrix1 &A, const Matrix2 &B,
				  MatrixCategories::RowColMatrixTag,
				  MatrixCategories::RowColMatrixTag) const
		{ return areEqualRow (A, B); }

	template <class Matrix> bool isZeroRow (const Matrix &v) const;
	template <class Matrix> bool isZeroCol (const Matrix &v) const;

	template <class Matrix>
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::RowMatrixTag) const
		{ return isZeroRow (A); }
	template <class Matrix>
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::ColMatrixTag) const
		{ return isZeroCol (A); }
	template <class Matrix>
	bool isZeroSpecialized (const Matrix &A, MatrixCategories::RowColMatrixTag) const
		{ return isZeroRow (A); }

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& addRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& addCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& addSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag) const
		{ return addRow (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& addSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag) const
		{ return addCol (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& addSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag) const
		{ return addRow (C, A, B); }

	template <class Matrix1, class Matrix2> Matrix1& addinRow (Matrix1 &A, const Matrix2 &B) const;
	template <class Matrix1, class Matrix2> Matrix1& addinCol (Matrix1 &A, const Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	inline Matrix1& addinSpecialized (Matrix1 &A, const Matrix2 &B,
					  MatrixCategories::RowMatrixTag,
					  MatrixCategories::RowMatrixTag) const
		{ return addinRow (A, B); }
	template <class Matrix1, class Matrix2>
	inline Matrix1& addinSpecialized (Matrix1 &A, const Matrix2 &B,
					  MatrixCategories::ColMatrixTag,
					  MatrixCategories::ColMatrixTag) const
		{ return addinCol (A, B); }
	template <class Matrix1, class Matrix2>
	inline Matrix1& addinSpecialized (Matrix1 &A, const Matrix2 &B,
					  MatrixCategories::RowColMatrixTag,
					  MatrixCategories::RowColMatrixTag) const
		{ return addinRow (A, B); }

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& subRow (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& subCol (Matrix1 &C, const Matrix2 &A, const Matrix3 &B) const;

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& subSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag) const
		{ return subRow (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& subSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag) const
		{ return subCol (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1& subSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag) const
		{ return subRow (C, A, B); }

	template <class Matrix1, class Matrix2> Matrix1& subinRow (Matrix1 &A, const Matrix2 &B) const;
	template <class Matrix1, class Matrix2> Matrix1& subinCol (Matrix1 &A, const Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	Matrix1& subinSpecialized (Matrix1 &A, const Matrix2 &B,
				   MatrixCategories::RowMatrixTag,
				   MatrixCategories::RowMatrixTag) const
		{ return subinRow (A, B); }
	template <class Matrix1, class Matrix2>
	Matrix1& subinSpecialized (Matrix1 &A, const Matrix2 &B,
				   MatrixCategories::ColMatrixTag,
				   MatrixCategories::ColMatrixTag) const
		{ return subinCol (A, B); }
	template <class Matrix1, class Matrix2>
	Matrix1& subinSpecialized (Matrix1 &A, const Matrix2 &B,
				   MatrixCategories::RowColMatrixTag,
				   MatrixCategories::RowColMatrixTag) const
		{ return subinRow (A, B); }

	template <class Matrix1, class Matrix2> Matrix1& negRow (Matrix1 &A, const Matrix2 &B) const;
	template <class Matrix1, class Matrix2> Matrix1& negCol (Matrix1 &A, const Matrix2 &B) const;

	template <class Matrix1, class Matrix2>
	inline Matrix1& negSpecialized (Matrix1 &A, const Matrix2 &B,
					MatrixCategories::RowMatrixTag,
					MatrixCategories::RowMatrixTag) const
		{ return negRow (A, B); }
	template <class Matrix1, class Matrix2>
	inline Matrix1& negSpecialized (Matrix1 &A, const Matrix2 &B,
					MatrixCategories::ColMatrixTag,
					MatrixCategories::ColMatrixTag) const
		{ return negCol (A, B); }
	template <class Matrix1, class Matrix2>
	inline Matrix1& negSpecialized (Matrix1 &A, const Matrix2 &B,
					MatrixCategories::RowColMatrixTag,
					MatrixCategories::RowColMatrixTag) const
		{ return negRow (A, B); }

	template <class Matrix> Matrix &neginRow (Matrix &A) const;
	template <class Matrix> Matrix &neginCol (Matrix &A) const;

	template <class Matrix>
	Matrix &neginSpecialized (Matrix &A, MatrixCategories::RowMatrixTag) const
		{ return neginRow (A); }
	template <class Matrix>
	Matrix &neginSpecialized (Matrix &A, MatrixCategories::ColMatrixTag) const
		{ return neginCol (A); }
	template <class Matrix>
	Matrix &neginSpecialized (Matrix &A, MatrixCategories::RowColMatrixTag) const
		{ return neginRow (A); }

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
				 MatrixCategories::ColMatrixTag) const
		{ return mulRowRowCol (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::ColMatrixTag) const
		{ return mulColRowCol (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::ColMatrixTag) const
		{ return mulRowRowCol (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag) const
		{ return mulRowRowRow (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag) const
		{ return mulColColCol (C, A, B); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &A, const Matrix3 &B,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag) const
		{ return mulRowRowCol (C, A, B); }

	template <class Matrix1, class Matrix2>
	Matrix1 &mulRow (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const;
	template <class Matrix1, class Matrix2>
	Matrix1 &mulCol (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a) const;

	template <class Matrix1, class Matrix2>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a,
				 MatrixCategories::RowMatrixTag,
				 MatrixCategories::RowMatrixTag) const
		{ return mulRow (C, B, a); }
	template <class Matrix1, class Matrix2>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a,
				 MatrixCategories::ColMatrixTag,
				 MatrixCategories::ColMatrixTag) const
		{ return mulCol (C, B, a); }
	template <class Matrix1, class Matrix2>
	Matrix1 &mulSpecialized (Matrix1 &C, const Matrix2 &B, const typename Field::Element &a,
				 MatrixCategories::RowColMatrixTag,
				 MatrixCategories::RowColMatrixTag) const
		{ return mulRow (C, B, a); }

	template <class Matrix> Matrix &mulinRow (Matrix &B, const typename Field::Element &a) const;
	template <class Matrix> Matrix &mulinCol (Matrix &B, const typename Field::Element &a) const;

	template <class Matrix>
	Matrix &mulinSpecialized (Matrix &B, const typename Field::Element &a,
				  MatrixCategories::RowMatrixTag) const
		{ return mulinRow (B, a); }
	template <class Matrix>
	Matrix &mulinSpecialized (Matrix &B, const typename Field::Element &a,
				  MatrixCategories::ColMatrixTag) const
		{ return mulinCol (B, a); }
	template <class Matrix>
	Matrix &mulinSpecialized (Matrix &B, const typename Field::Element &a,
				  MatrixCategories::RowColMatrixTag) const
		{ return mulinRow (B, a); }

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
		{ return axpyinRowRowCol (Y, A, X); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
				    MatrixCategories::ColMatrixTag,
				    MatrixCategories::RowMatrixTag,
				    MatrixCategories::ColMatrixTag) const
		{ return axpyinColRowCol (Y, A, X); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
				    MatrixCategories::RowColMatrixTag,
				    MatrixCategories::RowMatrixTag,
				    MatrixCategories::ColMatrixTag) const
		{ return axpyinRowRowCol (Y, A, X); }
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
				    MatrixCategories::RowMatrixTag,
				    MatrixCategories::RowMatrixTag,
				    MatrixCategories::RowMatrixTag) const
		{ return axpyinRowRowRow (Y, A, X); }
	
	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
				    MatrixCategories::ColMatrixTag,
				    MatrixCategories::ColMatrixTag,
				    MatrixCategories::ColMatrixTag) const
	{ return axpyinColColCol (Y, A, X); }

	template <class Matrix1, class Matrix2, class Matrix3>
	Matrix1 &axpyinSpecialized (Matrix1 &Y, const Matrix2 &A, const Matrix3 &X,
				    MatrixCategories::RowColMatrixTag,
				    MatrixCategories::RowColMatrixTag,
				    MatrixCategories::RowColMatrixTag) const
		{ return axpyinRowRowCol (Y, A, X); }

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::SparseSequenceVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::SparseAssociativeVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulRowSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::SparseParallelVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseVectorTag,
				    VectorCategories::DenseVectorTag) const
	{ return mulColDense (_VD, w, A, v); } 
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseVectorTag,
				    VectorCategories::SparseSequenceVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseVectorTag,
				    VectorCategories::SparseAssociativeVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				    VectorCategories::DenseVectorTag,
				    VectorCategories::SparseParallelVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	inline Vector1 &mulColSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
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

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 MatrixCategories::RowMatrixTag) const
		{ return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 MatrixCategories::ColMatrixTag) const
		{ return mulColSpecialized (w, A, v,
					    typename VectorTraits<Vector1>::VectorCategory (),
					    typename VectorTraits<Vector2>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &mulSpecialized (Vector1 &w, const Matrix &A, const Vector2 &v,
				 MatrixCategories::RowColMatrixTag) const
		{ return mulRowSpecialized (w, A, v, typename VectorTraits<Vector1>::VectorCategory ()); }

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::SparseSequenceVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::SparseAssociativeVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinRowSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::SparseParallelVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::DenseVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::SparseSequenceVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::SparseAssociativeVectorTag) const;
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinColSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				       VectorCategories::SparseParallelVectorTag) const;

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				    MatrixCategories::RowMatrixTag) const
		{ return axpyinRowSpecialized (y, A, x, typename VectorTraits<Vector1>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				    MatrixCategories::ColMatrixTag) const
		{ return axpyinColSpecialized (y, A, x, typename VectorTraits<Vector1>::VectorCategory ()); }
	template <class Vector1, class Matrix, class Vector2>
	Vector1 &axpyinSpecialized (Vector1 &y, const Matrix &A, const Vector2 &x,
				    MatrixCategories::RowColMatrixTag) const
		{ return axpyinRowSpecialized (y, A, x, typename VectorTraits<Vector1>::VectorCategory ()); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsByRow (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsByCol (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowColMatrixTag) const
		{ return permuteRowsByCol (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowMatrixTag) const
		{ return permuteRowsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteRowsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::ColMatrixTag) const
		{ return permuteRowsByCol (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsByRow (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsByCol (Matrix   &A,
					 Iterator  P_start,
					 Iterator  P_end) const;

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowColMatrixTag) const
		{ return permuteColsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::RowMatrixTag) const
		{ return permuteColsByRow (A, P_start, P_end); }

	template <class Matrix, class Iterator>
	inline Matrix &permuteColsSpecialized (Matrix   &A,
					       Iterator  P_start,
					       Iterator  P_end,
					       MatrixCategories::ColMatrixTag) const
		{ return permuteColsByCol (A, P_start, P_end); }

	Field         _F;
	VectorDomain<Field>  _VD;
};

}


#include "linbox/matrix/matrix-domain.inl"

#endif // __LINBOX_matrix_domain_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
