/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/direct-sum.h
 * David Saunders
 *
 * See COPYING about license
 */

#ifndef __DIRECT_SUM_H
#define __DIRECT_SUM_H

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/null-matrix.h"
#include "linbox/vector/vector-traits.h"
//#include "linbox/vector/subvector.h"

namespace LinBox
{

	/** @memo If C = DirectSum(A, B) and y = xA and z = wB, then (y,z) = (x,w)C.
	 * @doc
	 * And similarly for apply. 
	 */
	template <class Field, class Vector>
	class DirectSum : public BlackboxArchetype<Vector>
	{
	    public:

		typedef BlackboxArchetype<Vector> Blackbox;

		/** Constructor from two black box matrices.
		 * This becomes direct sum of A and B.
		 * They may be rectangular.  
		 * @param A, B:  black box matrices over a common field.
		 */
		DirectSum(const	Blackbox& A, const Blackbox& B)
			: _Ap(&A), _Bp(&B)
		{}

		/** Constructor from two black box matrix pointers.
		 * This becomes direct sum of A and B.
		 * They may be rectangular.  They must be over the same field (or ring). 
		 * @param A_ptr pointer to black box matrix A.
		 * @param B_ptr pointer to black box matrix B.
		 */
		DirectSum(): _Ap(& DirectSum<Field,Vector>::_NullMatrix),
			     _Bp(& DirectSum<Field,Vector>::_NullMatrix) 
		{}



		DirectSum(const	Blackbox* Ap, const Blackbox* Bp)
			: _Ap(Ap), _Bp(Bp)
		{}

		/// Copy constructor.
		DirectSum (const DirectSum<Field,Vector>& M) 
			: _Ap (M._Ap), _Bp (M._Bp)
		{}

		/// Destructor
		~DirectSum (void)
		{}

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Required by abstract base class.
		 */
		Blackbox* clone () const
		{ return new DirectSum(*this); }
		//	{ return new DirectSum (*this); }

		Vector& apply (Vector& y, const Vector& x) const
		{
		/* FIXME: I want to use subvectors to avoid copying and memory allocation, but problems with it...

			const Subvector<Vector> x1(x.begin(), x.begin() + _Ap->coldim());
			const Subvector<Vector> x2(x.begin() + _Ap->coldim(), x.end());
			Subvector<Vector> y1(y.begin(), y.begin() + _Ap->rowdim());
			Subvector<Vector> y2(y.begin() + _Ap->rowdim(), y.end());
		*/
			if (x.size() == 0) return y;  // Null matrix

			Vector xA(_Ap->coldim());
			Vector yA(_Ap->rowdim());
			std::copy (x.begin(), x.begin() + _Ap->coldim(), xA.begin());
			_Ap->apply (yA, xA);
			std::copy (yA.begin(), yA.end(), y.begin());

			Vector xB(_Bp->coldim());
			Vector yB(_Bp->rowdim());
			std::copy (x.begin() + _Ap->coldim(), x.end(), xB.begin());
			_Bp->apply (yB, xB);
			std::copy (yB.begin(), yB.end(), y.begin() + _Ap->rowdim());
			/*

			y.resize(_Ap->rowdim());
			x.resize(_Ap->coldim());
			_Ap->apply(y,x);
			*/
			return y;
		}

		Vector& applyTranspose (Vector& y, const Vector& x) const
		{
		/* FIXME: I want to use subvectors to avoid copying and memory allocation, but problems with it...
			const Subvector<Vector> x1(x.begin(), x.begin() + _Ap->rowdim());
			const Subvector<Vector> x2(x.begin() + _Ap->rowdim(), x.end());
			Subvector<Vector> y1(y.begin(), y.begin() + _Ap->coldim());
			Subvector<Vector> y2(y.begin() + _Ap->coldim(), y.end());

			Vector local_x(x1.size());
			Vector local_y(y1.size());
			copy (x1.begin(), x1.end(), local_x.begin());
			_Ap->applyTranspose (local_y, local_x);
			copy (local_y.begin(), local_y.end(), y1.begin());

			local_x.resize(x1.size());
			local_y.resize(y2.size());
			copy (x2.begin(), x2.end(), local_x.begin());
			_Bp->applyTranspose (local_y, local_x);
			copy (local_y.begin(), local_y.end(), y2.begin());

			//_Ap->applyTranspose (y1, x1);
			//_Bp->applyTranspose (y2, x2);

		*/
			if (x.size() == 0 ) return y;
			Vector xAT(_Ap->rowdim());
			Vector yAT(_Ap->coldim());
			std::copy (x.begin(), x.begin() + _Ap->rowdim(), xAT.begin());
			_Ap->apply (yAT, xAT);
			std::copy (yAT.begin(), yAT.end(), y.begin());

			Vector xBT(_Bp->rowdim());
			Vector yBT(_Bp->coldim());
			std::copy (x.begin() + _Ap->rowdim(), x.end(), xBT.begin());
			_Bp->apply (yBT, xBT);
			std::copy (yBT.begin(), yBT.end(), y.begin() + _Ap->rowdim());


			return y;
		}

		inline size_t rowdim (void) const
		{

			return _Ap->rowdim () + _Bp->rowdim ();

		}
    
		inline size_t coldim(void) const 
		{
			return _Ap->coldim () + _Bp->coldim ();
		}

	    protected:
		// the direct summands
		const Blackbox* _Ap;
		const Blackbox* _Bp; 
		static const NullMatrix<Vector>  _NullMatrix;
		

	}; // template <Vector> class DirectSum

	template<class Field, class Vector>	
	const NullMatrix<Vector> DirectSum<Field,Vector>::_NullMatrix;
	
}; // namespace LinBox

#endif // __DIRECT_SUM_H
