/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/compose.h
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

#include "linbox/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox compose matrix.
	 * This is a class that multiplies two matrices by implementing an 
	 * apply method that calls the apply methods of both of the consituent 
	 * matrices.
	 *
	 * This class, like the Black Box archetype from which it is derived, 
	 * is templatized by the vector type to which the matrix is applied.  
	 * Both constituent matrices must also use this same vector type.
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Vector>
	class Compose : public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;

		/** Constructor from two black box matrices.
		 * This constructor creates a matrix that is a product of two black box
		 * matrices: A*B.
		 * @param A_ptr pointer to black box matrix A.
		 * @param B_ptr pointer to black box matrix B.
		 */
		Compose (Blackbox *A_ptr, Blackbox *B_ptr)
		{
			linbox_check (A_ptr != (Blackbox *) 0);
			linbox_check (B_ptr != (Blackbox *) 0);
			linbox_check (A_ptr->coldim () == B_ptr->rowdim ());

			// create new copies of matrices in dynamic memory
			_A_ptr = A_ptr->clone ();
			_B_ptr = B_ptr->clone ();
			_z.resize (_B.ptr->rowdim ());
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Compose (const Compose<Vector>& M)
		{
			// create new copies of matrices in dynamic memory
			_A_ptr = M._A_ptr->clone ();
			_B_ptr = M._B_ptr->clone ();
			_z.resize (_B.ptr->rowdim ());
		}

		/// Destructor
		~Compose (void)
		{
			if (_A_ptr != (Blackbox *) 0) delete _A_ptr;
			if (_A_ptr != (Blackbox *) 0) delete _B_ptr;
		}

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the Blackbox_archetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox* clone () const
			{ return new Compose (*this); }

		/** Application of BlackBox matrix.
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

		/** Application of BlackBox matrix transpose.
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

		/** Retreive row dimensions of BlackBox matrix.
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
    
		/** Retreive column dimensions of BlackBox matrix.
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
		Blackbox *_A_ptr;
		Blackbox *_B_ptr;
		// local intermediate vector
		Vector _z;

	}; // template <Vector> class Compose

} // namespace LinBox

#endif // __COMPOSE
