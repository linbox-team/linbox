/* -*- mode: C++; style: linux -*- */

/* linbox/blackbox/sparse0.h    (Formerly sparse-matrix.h)
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 */

#ifndef __SPARSE0_H
#define __SPARSE0_H

#include "linbox/blackbox/sparse0-aux.h"
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
	 * It has two base classes.  One, BlackboxArchetype, is an abstract base
	 * class which ensures it adheres to the common object interface.  The
	 * other, sparsemat, contains the mathematical functions and objects
	 * necessary to implement the matrix.
	 *
	 * @param Field \Ref{LinBox} field
	 * @param Row	  \Ref{LinBox} sparse vector implementation for rows of matrix
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Field, class Row, class Vector>
	class SparseMatrix0
		: public BlackboxArchetype<Vector>, 
	public SparseMatrix0Aux<Field, Row, Vector>
	{
	public:
 
		/** Constructor from sparsemat_aux<Field, Row, Vector>.
		 * @param A constant reference to sparsemat object
		 */
		SparseMatrix0 (const SparseMatrix0Aux<Field, Row, Vector>& A) 
			: SparseMatrix0Aux<Field, Row, Vector> (A) {}

		/** Constructor.
		 * Note: the copy constructor and operator= will work as intended
		 *       because of STL's container design
		 * @param  F  the field of entries; passed so that a possible paramter 
		 *            such as a modulus is known to the matrix.
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		SparseMatrix0 (const Field& F, size_t m, size_t n)
			: SparseMatrix0Aux<Field, Row, Vector> (F, m, n) {}
    
		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		BlackboxArchetype<Vector>* clone () const 
			{ return new SparseMatrix0 (*this); }

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& apply (Vector& y, const Vector& x) const
			{ return SparseMatrix0Aux<Field, Row, Vector>::apply (y, x); }

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose (Vector& y, const Vector& x) const
			{ return SparseMatrix0Aux<Field, Row, Vector>::applyTranspose (y, x); }

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const 
			{ return SparseMatrix0Aux<Field, Row, Vector>::get_rowdim (); } 
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const 
			{ return SparseMatrix0Aux<Field, Row, Vector>::get_coldim (); } 

	}; // SparseMatrix0<Field> 

} // namespace LinBox

#endif // __SPARSE0_H
