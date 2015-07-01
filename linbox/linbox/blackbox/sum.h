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

#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"
#include <linbox/blackbox/blackbox-interface.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief blackbox of a matrix sum without copying.  

\ingroup blackbox
         * Adds only at apply time.
	 * Given two black boxes A and B of the same dimensions, form a black
	 * box representing A+B, i.e., Sum(A,B)x=(A+B)x=Ax+Bx
	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Blackbox1, class Blackbox2>
	class Sum : public BlackboxInterface 
	{
	    public:


		typedef typename Blackbox1::Field Field;
		typedef typename Blackbox1::Element Element;


		/** Constructor from black box matrices.
		 * This constructor creates a matrix that is the sum,
		 * A + B, of black box matrices.
		 * @param A, B:  black box matrices.
		 */
		Sum (const Blackbox1 &A, const Blackbox2 &B)
			: _A_ptr(&A), _B_ptr(&B)
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
		Sum (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr)
			: _A_ptr(A_ptr), _B_ptr(B_ptr)
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A_ptr != 0);
			linbox_check (B_ptr != 0);
			linbox_check (A_ptr->coldim () == B_ptr->coldim ());
			linbox_check (A_ptr->rowdim () == B_ptr->rowdim ());

			VectorWrapper::ensureDim (_z1, A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, A_ptr->coldim ());
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Sum (const Sum<Blackbox1, Blackbox2> &M)
			: _A_ptr (M._A_ptr), _B_ptr (M._B_ptr)
		{
			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/// Destructor
		~Sum (void)
		{
		}

		/** Application of BlackBox matrix.
		 * y= (A+B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		inline OutVector &apply (OutVector &y, const InVector &x) const
		{
			if (_A_ptr != 0 && _B_ptr != 0) {
				VectorDomain<Field> VD (field());
				typename Field::Element one;

				field().init (one, 1);
				_A_ptr->apply (y, x);
				_B_ptr->apply (_z1, x);
				VD.axpyin (y, one, _z1);
			}

			return y;
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A+B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		inline OutVector &applyTranspose (OutVector &y, const InVector &x) const
		{
			if (_A_ptr != 0 && _B_ptr != 0) {
				VectorDomain<Field> VD (field());
				typename Field::Element one;

				field().init (one, 1);
				_A_ptr->applyTranspose (y, x);
				_B_ptr->applyTranspose (_z2, x);
				VD.axpyin (y, one, _z2);
			}

			return y;
		}

            template<typename _Tp1, typename _Tp2 = _Tp1> 
            struct rebind                           
            { typedef Sum<typename Blackbox1::template rebind<_Tp1>::other, typename Blackbox2::template rebind<_Tp2>::other> other; };

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


		const Field& field() const { return _A_ptr -> field(); }

	    protected:

		// use a copy of the input field for faster performance (no pointer dereference).

		const Blackbox1       *_A_ptr;
		const Blackbox2       *_B_ptr;

		mutable std::vector<Element>  _z1;
		mutable std::vector<Element>  _z2;

	}; // template <Field, Vector> class Sum

} // namespace LinBox

#endif // __SUM_H
