/* linbox/blackbox/matrix-blackbox.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * --------------------------------------------------------
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

#ifndef __LINBOX_matrix_blackbox_H
#define __LINBOX_matrix_blackbox_H

#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/field/rebind.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief Matrix black box
	 *
	 \ingroup blackbox
	 * This class wraps a matrix meeting the @ref MatrixArchetype  interface into a
	 * black box meeting the @ref BlackboxArchetype interface. It uses
	 * @ref MatrixDomain to implement \c apply and \c applyTranspose.
	 */

	template <class _Field, class _Matrix, class _Vector = typename LinBox::Vector<_Field>::Dense>
	class MatrixBlackbox : public BlackboxArchetype {
	public:

		typedef _Vector                 Vector;
		typedef _Matrix                 Matrix;
		typedef _Field                  Field;
		typedef typename Field::Element Element;

		/** Constructor.
		 *
		 * Builds a black box for the matrix given by \p rep over the field
		 * \p F
		 *
		 * @param  F  Field over which entries exist
		 * @param  rep  Matrix from which to construct the black box
		 */
		MatrixBlackbox (const Field &F, Matrix &Rep) :
			_field (&F), _MD (F), _matA (Rep)
		{}

		/** Constructor with size
		 *
		 * Builds a black box with the given dimensions
		 *
		 * @param  F  Field over which entries exist
		 * @param  m  Row dimension
		 * @param  n  Column dimension
		 */
		MatrixBlackbox (const Field &F, size_t m, size_t n) :
			_field (&F), _MD (F), _matA (F, m, n)
		{}

		/** Constructor
		 *
		 * Builds a black box, using the vector stream \p stream to fill in
		 * its entries
		 *
		 * @param  F  Field over which entries exist
		 * @param  stream  Stream with which to generate row vectors
		 */
		template <class Row>
		MatrixBlackbox (const Field &F, VectorStream<Row> &stream) :
			_field (&F), _MD (F), _matA (stream)
		{}

		/** Copy constructor
		*/
		MatrixBlackbox (const MatrixBlackbox &B) :
			_field (B._field), _MD (*(B._field)), _matA (B._matA)
		{}

		/** Destructor. */
		~MatrixBlackbox () {}


		template<typename _Tp1,
		typename _Mat1 = typename Matrix::template rebind<_Tp1>::other,
		typename _Vect1 = typename Rebind<Vector, _Tp1>::other >
		struct rebind {
			typedef MatrixBlackbox<_Tp1, _Mat1, _Vect1> other;
		};



		/** Create a clone of the matrix
		*/
		inline BlackboxArchetype *clone () const
		{ return new MatrixBlackbox (*this); }

		/** Generic matrix-vector product
		 * \f$ y = A x\f$.
		 * @return reference to output vector y
		 * @param  x input vector
		 * @param y
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &apply (Vector1 &y, const Vector2 &x) const
		{ return _MD.vectorMul (y, _matA, x); }

		/** Matrix-vector product
		 * \f$y = A x\f$.
		 * @return reference to output vector y
		 * @param  x input vector
		 * @param y
		 */
		inline Vector &apply (Vector &y, const Vector &x) const
		{ return apply<Vector, Vector> (y, x); }

		/** Generic transpose matrix-vector product
		 * \f$y = A^T x\f$.
		 * @return reference to output vector y
		 * @param  x input vector
		 * @param y
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const
		{ return _MD.vectorMul (y, TransposeMatrix<const Matrix> (_matA), x); }

		/** Transpose matrix-vector product
		 * \f$y = A^T x\f$.
		 * @return reference to output vector y
		 * @param  x input vector
		 * @param y
		 */
		inline Vector &applyTranspose (Vector &y, const Vector &x) const
		{ return applyTranspose<Vector, Vector> (y, x); }

		/** Retreive row dimensions of Sparsemat matrix.
		 * @return integer number of rows of SparseMatrix0Base matrix.
		 */
		inline size_t rowdim () const
		{ return _matA.rowdim (); }

		/** Retreive column dimensions of Sparsemat matrix.
		 * @return integer number of columns of SparseMatrix0Base matrix.
		 */
		inline size_t coldim () const
		{ return _matA.coldim (); }

		/** Read the matrix from a stream
		 * @param is Input stream from which to read the matrix
		 * @return Reference to input stream
		 */
		inline std::istream &read (std::istream &is)
		{ return _MD.read (is, _matA); }

		/** Write the matrix to a stream
		 * @param os Output stream to which to write the matrix
		 * @return Reference to output stream
		 */
		inline std::ostream &write (std::ostream &os) const
		{ return _MD.write (os, _matA); }

		/** Return a reference to the base field
		*/
		inline const Field &field () const { return *_field;}

		/** Return a reference to the underlying representation
		*/
		inline Matrix &rep () { return _matA; }

	private:

		const Field         *_field;      // Field used for all arithmetic
		MatrixDomain<Field>  _MD;     // Matrix domain for matrix-vector
		// operations
		Matrix               _matA;      // Underlying matrix representation
	};

} // namespace LinBox

#endif // __LINBOX_matrix_blackbox_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

