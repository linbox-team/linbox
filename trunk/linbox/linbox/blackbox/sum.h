/* -*- mode: c; style: linux -*- */

/* linbox/blackbox/sum.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __TRANSPOSE_H
#define __TRANSPOSE_H

#include "linbox/blackbox/archetype.h"
#include "linbox/field/vector-domain.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Given two black boxes A and B of the same dimensions, form a black
	 * box representing A+B, i.e., Sum(A,B)x=(A+B)x=Ax+Bx
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Field, class Vector>
	class Sum : public BlackboxArchetype<Vector>
	{
	    public:

		typedef BlackboxArchetype<Vector> Blackbox;

		/** Constructor from a black box.
		 * This constructor creates a matrix that the transpose of a black box
		 * matrix A
		 * @param A_ptr pointer to black box matrix.
		 */
		Sum (const Field &F, const Blackbox *A_ptr, const Blackbox *B_ptr)
			: _F (F)
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A_ptr != 0);
			linbox_check (B_ptr != 0);
			linbox_check (A_ptr->coldim () == B_ptr->coldim ());
			linbox_check (A_ptr->rowdim () == B_ptr->rowdim ());

			_A_ptr = A_ptr->clone ();
			_B_ptr = B_ptr->clone ();
			VectorWrapper::ensureDim (_z1, A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, A_ptr->coldim ());
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Sum (const Sum<Field, Vector> &M)
			: _F (M._F), _A_ptr (M._A_ptr->clone ()), _B_ptr (M._B_ptr->clone ())
		{
			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/// Destructor
		~Sum (void)
		{
			delete _A_ptr;
			delete _B_ptr;
		}

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the BlackboxArchetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox *clone () const
			{ return new Sum (*this); }

		/** Application of BlackBox matrix.
		 * y= (A+B)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		inline Vector &apply (Vector &y, const Vector &x) const
		{
			if (_A_ptr != 0 && _B_ptr != 0) {
				VectorDomain<Field> VD (_F);
				typename Field::Element one;

				_F.init (one, 1);
				_A_ptr->apply (y, x);
				_B_ptr->apply (_z1, x);
				VD.axpyin (y, one, _z1);
			}

			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A+B)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		inline Vector &applyTranspose (Vector &y, const Vector &x) const
		{
			if (_A_ptr != 0 && _B_ptr != 0) {
				VectorDomain<Field> VD (_F);
				typename Field::Element one;

				_F.init (one, 1);
				_A_ptr->applyTranspose (y, x);
				_B_ptr->applyTranspose (_z2, x);
				VD.axpyin (y, one, _z2);
			}

			return y;
		}

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
			{ return _A_ptr->rowdim (); }
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const 
			{ return _A_ptr->coldim (); }

	    private:

		const Field    &_F;

		Blackbox       *_A_ptr;
		Blackbox       *_B_ptr;

		mutable Vector  _z1;
		mutable Vector  _z2;

	}; // template <Field, Vector> class Sum

} // namespace LinBox

#endif // __TRANSPOSE_H
