/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/sparse-matrix.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __SPARSE_MATRIX_H
#define __SPARSE_MATRIX_H

#include "linbox/blackbox/sparse-matrix-aux.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox sparse matrix template.
	 * This is a class of sparse matrices templatized by the 
	 * {@link Fields field} in
	 * which the elements reside.  The matrix itself is stored as an
	 * STL vector of \Ref{LinBox} sparse vectors of integers and field elements.
	 * Each sparse vector corresponds to one row of the matrix, 
	 * and each pair (j, a) in sparse vector i corresponds to the (i,j) 
	 * entry of the matrix.
	 *
	 * The class is conforms to the 
	 * {@link Archetypes archetype} for \Ref{Blackbox Matrices}.
	 *
	 * It has two base classes.  One, Blackbox_archetype, is an abstract base
	 * class which ensures it adheres to the common object interface.  The
	 * other, sparsemat, contains the mathematical functions and objects
	 * necessary to implement the matrix.
	 *
	 * @param Field \Ref{LinBox} field
	 * @param Row	  \Ref{LinBox} sparse vector implementation for rows of matrix
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Field, class Row, class Vector>
	class SparseMatrix
		: public Blackbox_archetype<Vector>, public SparseMatrixAux<Field, Row, Vector>
	{
	public:
 
		/** Constructor from sparsemat_aux<Field, Row, Vector>.
		 * @param A constant reference to sparsemat object
		 */
		SparseMatrix (const SparseMatrixAux<Field, Row, Vector>& A) 
			: SparseMatrixAux<Field, Row, Vector> (A) {}

		/** Constructor.
		 * Note: the copy constructor and operator= will work as intended
		 *       because of STL's container design
		 * @param  F  the field of entries; passed so that a possible paramter 
		 *            such as a modulus is known to the matrix.
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		SparseMatrix (const Field& F, size_t m, size_t n)
			: SparseMatrixAux<Field, Row, Vector> (F, m, n) {}
    
		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the Blackbox_archetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox_archetype<Vector>* clone () const 
			{ return new SparseMatrix (*this); }

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& apply (Vector& y, const Vector& x) const
			{ return SparseMatrixAux<Field, Row, Vector>::apply (y, x); }

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose (Vector& y, const Vector& x) const
			{ return SparseMatrixAux<Field, Row, Vector>::applyTranspose (y, x); }

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const 
			{ return SparseMatrixAux<Field, Row, Vector>::get_rowdim (); } 
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const 
			{ return SparseMatrixAux<Field, Row, Vector>::get_coldim (); } 

	}; // SparseMatrix<Field> 

} // namespace LinBox

#endif // __SPARSE_MATRIX_H
