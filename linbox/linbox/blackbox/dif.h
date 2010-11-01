/* linbox/blackbox/dif.h
 * Copyright (c) 2010 LinBox
 * transmuted from linbox/blackbox/sum.h by bds
 * Modified by JG Dumas
 *
 * It will be desirable to keep sum.h and dif.h in sync.
 * 
 * Time-stamp: <18 Jun 10 15:28:48 Jean-Guillaume.Dumas@imag.fr> 
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_dif_H
#define __LINBOX_dif_H

#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"
#include <linbox/blackbox/blackbox-interface.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief Blackbox of a difference: C := A - B, i.e. Cx = Ax - Bx.

	 * {\bf Template parameters:} 
	 *     Field is the class of the entry domain, 
	 *     Vector is a LinBox dense or sparse vector of field elements class.
\ingroup blackbox
	 */
	template <class _Blackbox1, class _Blackbox2>
	class Dif : public BlackboxInterface 
	{
                typedef Dif<_Blackbox1, _Blackbox2> Self_t;
	    public:
		typedef _Blackbox1 Blackbox1;
		typedef _Blackbox2 Blackbox2;
	       
		typedef typename Blackbox1::Field Field;
		typedef typename Blackbox1::Element Element;

		/** Build this as A - B from blackboxes A, B.
		 * A and B must have the same shape and be over the same field.
		 * Their data is not copied.  A subsequent change to one of them also changes
		 * this difference.
		 */
		Dif (const Blackbox1 &A, const Blackbox2 &B) : _A_ptr(&A), _B_ptr(&B)		     
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A.coldim () == B.coldim ());
			linbox_check (A.rowdim () == B.rowdim ());

			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/** Build this as A - B from blackbox pointers A_ptr, B_ptr.
		 * The two matrices must have the same shape and be over the same field.
		 * Their data {\it is} copied.  I don't know why.
		 */
		Dif (const Blackbox1 *A_ptr, const Blackbox2 *B_ptr) : _A_ptr(A_ptr),_B_ptr(B_ptr)	       
		{
			// create new copies of matrices in dynamic memory
			linbox_check (A_ptr != 0);
			linbox_check (B_ptr != 0);
			linbox_check (A_ptr->coldim () == B_ptr->coldim ());
			linbox_check (A_ptr->rowdim () == B_ptr->rowdim ());

// 			_A_ptr = A_ptr->clone ();
// 			_B_ptr = B_ptr->clone ();
			VectorWrapper::ensureDim (_z1, A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, A_ptr->coldim ());
		}

		/** Makes a deep copy.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Dif (const Dif<Blackbox1, Blackbox2> &M)
			: _A_ptr (M._A_ptr), _B_ptr (M._B_ptr)
		{
			VectorWrapper::ensureDim (_z1, _A_ptr->rowdim ());
			VectorWrapper::ensureDim (_z2, _A_ptr->coldim ());
		}

		/// Destructor
		~Dif (void)
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

			if ((_A_ptr != 0) && (_B_ptr != 0)) {

				VectorDomain<Field> VD (_A_ptr->field());

				_A_ptr->apply (y, x);
				_B_ptr->apply (_z1, x);

				VD.subin(y, _z1);
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
			if ((_A_ptr != 0) && (_B_ptr != 0)) {
				VectorDomain<Field> VD (_A_ptr->field());
				_A_ptr->applyTranspose (y, x);
				_B_ptr->applyTranspose (_z2, x);
				VD.subin (y, _z2);
			}

			return y;
		}

            template<typename _Tp1, typename _Tp2 = _Tp1>
            struct rebind
            { 
                typedef Dif<
                    typename Blackbox1::template rebind<_Tp1>::other,
                    typename Blackbox2::template rebind<_Tp2>::other
                > other;
                
    		void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                    typename Blackbox1::template rebind<_Tp1> () ( *(Ap._A_ptr), *(A._A_ptr), F);
                    typename Blackbox2::template rebind<_Tp2> () ( *(Ap._B_ptr), *(A._B_ptr), F);
                }

            };



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

		const Field& field() const {return _A_ptr->field();}

	    private:

		const Blackbox1       *_A_ptr;
		const Blackbox2       *_B_ptr;

		mutable std::vector<Element>  _z1;
		mutable std::vector<Element>  _z2;

	}; 

} // namespace LinBox

#endif // __LINBOX_dif_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
