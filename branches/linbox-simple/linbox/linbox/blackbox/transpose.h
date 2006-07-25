/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __TRANSPOSE_H
#define __TRANSPOSE_H

#include <linbox/blackbox/blackbox-interface.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** \brief transpose matrix without copying.

\ingroup blackbox

	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Blackbox>
	class Transpose : public BlackboxInterface
	{

	    public:
		typedef Blackbox Blackbox_t;
		typedef Transpose<Blackbox> Self_t;

		typedef typename Blackbox::Field Field;
		typedef typename Blackbox::Element Element;

		/** Constructor from a black box.
		 * This constructor creates a matrix that the transpose of a black box
		 * matrix A
		 * @param A_ptr pointer to black box matrix.
		 */
		Transpose (const Blackbox& A) : _A_ptr(&A){}

		Transpose (const Blackbox *A_ptr): _A_ptr(A_ptr)
		{
			// create new copies of matrices in dynamic memory
			//linbox_check (A_ptr != NULL);
			//_A_ptr = A_ptr->clone ();
		}

		/** Copy constructor.
		 * Creates new black box objects in dynamic memory.
		 * @param M constant reference to compose black box matrix
		 */
		Transpose (const Transpose<Blackbox> &M) : _A_ptr(M._A_ptr)
		{
			// create new copies of matrices in dynamic memory
			//linbox_check (M._A_ptr != NULL);
			//_A_ptr = M._A_ptr->clone ();
		}

		/// Destructor
		~Transpose (void)
		{
		}

                    template<typename _Tp1> 
                    struct rebind 
                    { 
                        typedef Transpose<typename Blackbox::template rebind<_Tp1>::other> other; 
                        void operator() (other *& Ap, const Self_t& A, const _Tp1& F) {
                            typename other::Blackbox_t * A1;
                            typename Blackbox_t::template rebind<_Tp1> () ( A1, *(A._A_ptr), F);
                            Ap = new other(A1);
                        }
                    };


		/** Application of BlackBox matrix.
		 * y= (A*B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &apply (Vector1 &y, const Vector2 &x) const
		{
			if (_A_ptr != 0) _A_ptr->applyTranspose (y, x);
			return y;
		}


		/** Application of BlackBox matrix transpose.
		 * y= transpose(A*B)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template <class Vector1, class Vector2>
		inline Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const
		{
			if (_A_ptr != 0) _A_ptr->apply (y, x);
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
				return _A_ptr->coldim ();
			else 
				return 0;
		}
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const 
		{
			if (_A_ptr != 0) 
				return _A_ptr->rowdim ();
			else 
				return 0;
		}
	       

		const Field& field() const {return _A_ptr->field();}
	    private:

		// Pointers to A and B matrices
		const Blackbox *_A_ptr;

	}; // template <Vector> class Transpose

} // namespace LinBox

#endif // __TRANSPOSE_H
