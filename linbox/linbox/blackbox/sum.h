/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sum.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __SUM_H
#define __SUM_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** @memo Blackbox of matrix sum without copying.  Add only at apply time.
	 * @doc Given two black boxes A and B of the same dimensions, form a black
	 * box representing A+B, i.e., Sum(A,B)x=(A+B)x=Ax+Bx
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Field, class Vector>
	class Sum : public BlackboxArchetype<Vector>
	{
	    public:

		typedef BlackboxArchetype<Vector> Blackbox;

		/** Constructor from black box matrices.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A, B:  black box matrices.
		 */
		Sum (const Field &F, const Blackbox &A, const Blackbox &B)
			: _F (F), _A_ptr(&A), _B_ptr(&B)
		{
			linbox_check (A.coldim () == B.coldim ());
			linbox_check (A.rowdim () == B.rowdim ());

			VectorWrapper::ensureDim (_z1, A.rowdim ());
			VectorWrapper::ensureDim (_z2, A.coldim ());
		}

		/** Constructor from black box pointers.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A_ptr, B_ptr:  pointers to black box matrices.
		 */
		Sum (const Field &F, const Blackbox *A_ptr, const Blackbox *B_ptr)
			: _F (F), _A_ptr(A_ptr), _B_ptr(B_ptr)
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A_ptr != 0);
			linbox_check (B_ptr != 0);
			linbox_check (A_ptr->coldim () == B_ptr->coldim ());
			linbox_check (A_ptr->rowdim () == B_ptr->rowdim ());

			// we won't copy the matrices
			//_A_ptr = A_ptr->clone ();
			//_B_ptr = B_ptr->clone ();
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
			//std::cerr << "Sum destructor called" << std::endl;
			//delete _A_ptr;
			//delete _B_ptr;
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

	    protected:

		// use a copy of the input field for faster performance (no pointer dereference).
		const Field    _F;
		//const Field    &_F;

		const Blackbox       *_A_ptr;
		const Blackbox       *_B_ptr;

		mutable Vector  _z1;
		mutable Vector  _z2;

	}; // template <Field, Vector> class Sum

} // namespace LinBox

#endif // __SUM_H
