/* linbox/blackbox/matrix-blackbox.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * --------------------------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_blackbox_H
#define __LINBOX_matrix_blackbox_H

#include <linbox/blackbox/blackbox-interface.h>
#include "linbox/blackbox/archetype.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/sparse.h"
#include "linbox/matrix/dense.h"
#include "linbox/field/rebind.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

/** \brief Matrix black box
 *
\ingroup blackbox
 * This class wraps a matrix meeting the @ref{MatrixArchetype} interface into a
 * black box meeting the @ref{BlackboxArchetype} interface. It uses
 * @ref{MatrixDomain} to implement @code{apply} and @code{applyTranspose}.
 */

template <class _Field, class _Matrix, class _Vector = typename LinBox::Vector<_Field>::Dense>
class MatrixBlackbox : public BlackboxArchetype
{
    public:

	typedef _Vector                 Vector;
	typedef _Matrix                 Matrix;
	typedef _Field                  Field;
	typedef typename Field::Element Element;

	/** Constructor.
	 *
	 * Builds a black box for the matrix given by @code{rep} over the field
	 * @code{F} 
	 *
	 * @param  F  Field over which entries exist
	 * @param  rep  Matrix from which to construct the black box
	 */
	MatrixBlackbox (const Field &F, Matrix &rep)
		: _F (F), _MD (F), _A (rep) {}

	/** Constructor with size
	 *
	 * Builds a black box with the given dimensions
	 *
	 * @param  F  Field over which entries exist
	 * @param  m  Row dimension
	 * @param  n  Column dimension
	 */
	MatrixBlackbox (const Field &F, size_t m, size_t n)
		: _F (F), _MD (F), _A (m, n) {}

	/** Constructor
	 *
	 * Builds a black box, using the vector stream @code{stream} to fill in
	 * its entries
	 *
	 * @param  F  Field over which entries exist
	 * @param  stream  Stream with which to generate row vectors
	 */
	template <class Row>
	MatrixBlackbox (const Field &F, VectorStream<Row> &stream)
		: _F (F), _MD (F), _A (stream) {}

	/** Copy constructor
	 */
	MatrixBlackbox (const MatrixBlackbox &B)
		: _F (B._F), _MD (B._F), _A (B._A) {}

	/** Destructor. */
	~MatrixBlackbox () {}


    template<typename _Tp1, 
             typename _Mat1 = typename Matrix::template rebind<_Tp1>::other, 
             typename _Vect1 = typename Rebind<Vector, _Tp1>::other >
    struct rebind
    { typedef MatrixBlackbox<_Tp1, _Mat1, _Vect1> other; };



	/** Create a clone of the matrix
	 */
	inline BlackboxArchetype *clone () const
		{ return new MatrixBlackbox (*this); }

	/** Generic matrix-vector product
	 * y = A x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	template <class Vector1, class Vector2>
	inline Vector1 &apply (Vector1 &y, const Vector2 &x) const
		{ return _MD.vectorMul (y, _A, x); }

	/** Matrix-vector product
	 * y = A x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	inline Vector &apply (Vector &y, const Vector &x) const
		{ return apply<Vector, Vector> (y, x); }

	/** Generic transpose matrix-vector product
	 * y = A^T x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	template <class Vector1, class Vector2>
	inline Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const
		{ return _MD.vectorMul (y, TransposeMatrix<const Matrix> (_A), x); }

	/** Transpose matrix-vector product
	 * y = A^T x.
	 * @return reference to output vector y
	 * @param  x input vector
	 */
	inline Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTranspose<Vector, Vector> (y, x); }

	/** Retreive row dimensions of Sparsemat matrix.
	 * @return integer number of rows of SparseMatrix0Base matrix.
	 */
	inline size_t rowdim () const
		{ return _A.rowdim (); }

	/** Retreive column dimensions of Sparsemat matrix.
	 * @return integer number of columns of SparseMatrix0Base matrix.
	 */
	inline size_t coldim () const
		{ return _A.coldim (); }

	/** Read the matrix from a stream
	 * @param is Input stream from which to read the matrix
	 * @return Reference to input stream
	 */
	inline std::istream &read (std::istream &is)
		{ return _MD.read (is, _A); }

	/** Write the matrix to a stream
	 * @param os Output stream to which to write the matrix
	 * @return Reference to output stream
	 */
	inline std::ostream &write (std::ostream &os) const
		{ return _MD.write (os, _A); }

	/** Return a reference to the base field
	 */
	inline const Field &field () const { return _F;}

	/** Return a reference to the underlying representation
	 */
	inline Matrix &rep () { return _A; }

    private:

	const Field         &_F;      // Field used for all arithmetic
	MatrixDomain<Field>  _MD;     // Matrix domain for matrix-vector
				      // operations
	Matrix               _A;      // Underlying matrix representation
};

} // namespace LinBox

#endif // __LINBOX_matrix_blackbox_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
