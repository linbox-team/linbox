/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/compose.h
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

#ifndef __COMPOSE_H
#define __COMPOSE_H

#include "linbox/blackbox/archetype.h"

#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** @memo Compose two blackboxes: C := AB, i.e. Cx := A(Bx).
	 * @doc
	 * This is a class that multiplies two matrices by implementing an 
	 * apply method that calls the apply methods of both of the consituent 
	 * matrices, one after the other.
	 *
	 * This class, like the Black Box archetype from which it is derived, 
	 * is templatized by the vector type to which the matrix is applied.  
	 * Both constituent matrices must also use this same vector type.
	 * For specification of the blackbox members see \Ref{BlackboxArchetype}.
	 * 
	 * {\bf Template parameter:} must meet the \Ref{Vector} requirement.
	 */
	template <class _Vector>
	class Compose : public BlackboxArchetype<_Vector>
	{
	    public:

		typedef _Vector Vector;
		typedef BlackboxArchetype<Vector> Blackbox;

		/** Constructor of C := A*B from blackbox matrices A and B.
		 * Build the product A*B of any two black box matrices of compatible dimensions.
		 * Requires A.coldim() equals B.rowdim().
		 */
		Compose (const Blackbox &A, const Blackbox &B)
			: _A_ptr(&A), _B_ptr(&B) 
		{
			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
		}

		/** Constructor of C := (*A_ptr)*(*B_ptr).
		 * This constructor creates a matrix that is a product of two black box
		 * matrices: A*B from pointers to them.
		 */
		Compose (const Blackbox *A_ptr, const Blackbox *B_ptr)
			: _A_ptr(A_ptr), _B_ptr(B_ptr)
		{
			linbox_check (A_ptr != (Blackbox *) 0);
			linbox_check (B_ptr != (Blackbox *) 0);
			linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

			// create new copies of matrices in dynamic memory
			//_A_ptr = A_ptr->clone ();
			//_B_ptr = B_ptr->clone ();

			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
		}

		/** Copy constructor.
		 * Copies the composed matrix (a small handle).  The underlying two matrices
		 * are not copied.
		 */
		Compose (const Compose<Vector>& M) 
			:_A_ptr ( M._A_ptr), _B_ptr ( M._B_ptr)
		{
			// create new copies of matrices in dynamic memory
			//_A_ptr = M._A_ptr->clone ();
			//_B_ptr = M._B_ptr->clone ();

			VectorWrapper::ensureDim (_z, _A_ptr->coldim ());
		}

		/// Destroy composition object, but not the underlying two matrices.
		virtual ~Compose (void)
		{
			//if (_A_ptr != (Blackbox *) 0) delete _A_ptr;
			//if (_A_ptr != (Blackbox *) 0) delete _B_ptr;
		}

		/*- Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox* clone () const
			{ return new Compose (*this); }

		/*- Application of BlackBox matrix.
		 * y= (A*B)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		inline Vector& apply (Vector& y, const Vector& x) const
		{
			if ((_A_ptr != 0) && (_B_ptr != 0)) {
				_B_ptr->apply (_z, x);
				_A_ptr->apply (y, _z);
			}

			return y;
		}

		/*- Application of BlackBox matrix transpose.
		 * y= transpose(A*B)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		inline Vector& applyTranspose (Vector& y, const Vector& x) const
		{
			if ((_A_ptr != 0) && (_B_ptr != 0)) {
				_A_ptr->applyTranspose (_z, x);
				_B_ptr->applyTranspose (y, _z);
			}

			return y;
		}

		/*- Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			if (_A_ptr != 0) 
				return _A_ptr->rowdim ();
			else 
				return 0;
		}
    
		/*- Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const 
		{
			if (_B_ptr != 0) 
				return _B_ptr->coldim ();
			else 
				return 0;
		}

	    private:

		// Pointers to A and B matrices
		const Blackbox *_A_ptr;
		const Blackbox *_B_ptr;

		// local intermediate vector
		mutable Vector _z;

	}; // template <Vector> class Compose

} // namespace LinBox

#endif // __COMPOSE
