/* linbox/blackbox/transpose.h
 * Copyright (C) 2001 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __LINBOX_squarize_H
#define __LINBOX_squarize_H

#include <linbox/blackbox/blackbox-interface.h>

#ifndef GIVMAX
#define GIVMAX(a,b) ((b)>(a)?(b):(a))
#endif

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief transpose matrix without copying.

\ingroup blackbox

	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Blackbox>
	class Squarize : public BlackboxInterface
	{

	    public:
		typedef Blackbox Blackbox_t;
		typedef Squarize<Blackbox> Self_t;

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;

		/** Constructor from a black box.
		 * This constructor creates a matrix that the transpose of a black box
		 * matrix A
		 * @param A_ptr pointer to black box matrix.
		 */
		Squarize (const Blackbox& A) : _A_ptr(&A){
			_A_ptr->field().init(_Zero,0UL);
		}

		Squarize (const Blackbox *A_ptr): _A_ptr(A_ptr)
		{
			_A_ptr->field().init(_Zero,0UL);
			// create new copies of matrices in dynamic memory
			//linbox_check (A_ptr != NULL);
			//_A_ptr = A_ptr->clone ();
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Squarize (const Squarize<Blackbox> &M) : _A_ptr(M._A_ptr)
		{
			_A_ptr->field().init(_Zero,0UL);
			// create new copies of matrices in dynamic memory
			//linbox_check (M._A_ptr != NULL);
			//_A_ptr = M._A_ptr->clone ();
		}

		/// Destructor
		~Squarize (void)
		{
		}

		template<typename _Tp1> 
		struct rebind 
		{ 
                        typedef Squarize<typename Blackbox::template rebind<_Tp1>::other> other; 
                        void operator() (other & Ap, const Self_t& A, const _Tp1& F) {
                                typename Blackbox_t::template rebind<_Tp1> () ( *(Ap._A_ptr), *(A._A_ptr), F);
                        }
		};

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &apply (Vector1 &y, const Vector2 &x) const
		{
			if (_A_ptr != 0) _A_ptr->apply (y, x);
			if (_A_ptr->rowdim () < y.size()) {
				for(typename Vector1::iterator yit=y.begin()+_A_ptr->rowdim ();
				    yit != y.end(); ++yit)
					*yit = _Zero;
			}
			return y;
		}


		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const
		{
			if (_A_ptr != 0) _A_ptr->applyTranspose (y, x);
			if (_A_ptr->coldim () < y.size()) {
				for(typename Vector1::iterator yit=y.begin()+_A_ptr->rcoldim ();
				    yit != y.end(); ++yit)
					*yit = _Zero;
}			return y;
		}

	protected:
		size_t maxsize() const {
			if (_A_ptr != 0) 
				return GIVMAX( _A_ptr->rowdim (), _A_ptr->coldim () );
			else 
				return 0;
		}		
	public:
		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const
		{
			return maxsize();
		}
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const 
		{
			return maxsize();
		}
	       


		const Field& field() const {return _A_ptr->field();}
	    protected:
		// Pointer to A matrix
		const Blackbox *_A_ptr;
		typename Field::Element _Zero;

	}; // template <Vector> class Squarize

} // namespace LinBox

#endif // __LINBOX_squarize_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
