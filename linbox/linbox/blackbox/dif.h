/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/dif.h
 * transmuted from linbox/blackbox/sum.h by bds
 *
 * It will be desirable to keep sum.h and dif.h in sync.
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __DIF_H
#define __DIF_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** @memo Blackbox of a difference: C := A - B, i.e. Cx = Ax - Bx.
	 * @doc 
	 * {\bf Template parameters:} 
	 *     Field is the class of the entry domain, 
	 *     Vector is a LinBox dense or sparse vector of field elements class.
	 */
	template <class Field, class Vector>
	class Dif : public BlackboxArchetype<Vector>
	{
	    public:

		typedef BlackboxArchetype<Vector> Blackbox;

		/** Build this as A - B from blackboxes A, B.
		 * A and B must have the same shape and be over the same field.
		 * Their data is not copied.  A subsequent change to one of them also changes
		 * this difference.
		 */
		Dif (const Field &F, const Blackbox &A, const Blackbox &B)
			: _F (F)
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A.coldim () == B.coldim ());
			linbox_check (A.rowdim () == B.rowdim ());

			_A_ptr = &A;
			_B_ptr = &B;
			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/** Build this as A - B from blackbox pointers A_ptr, B_ptr.
		 * The two matrices must have the same shape and be over the same field.
		 * Their data {\it is} copied.  I don't know why.
		 */
		Dif (const Field &F, const Blackbox *A_ptr, const Blackbox *B_ptr)
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

		/** Makes a deep copy.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Dif (const Dif<Field, Vector> &M)
			: _F (M._F), _A_ptr (M._A_ptr->clone ()), _B_ptr (M._B_ptr->clone ())
		{
			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/// Destructor
		~Dif (void)
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
			{ return new Dif (*this); }

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
				_A_ptr->apply (y, x);
				_B_ptr->apply (_z1, x);
				VD.subin(y, _z1);
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
				_A_ptr->applyTranspose (y, x);
				_B_ptr->applyTranspose (_z2, x);
				VD.subin (y, _z2);
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

		const Field    _F;

		Blackbox       *_A_ptr;
		Blackbox       *_B_ptr;

		mutable Vector  _z1;
		mutable Vector  _z2;

	}; // template <Field, Vector> class Dif

} // namespace LinBox

#endif // __DIF_H
