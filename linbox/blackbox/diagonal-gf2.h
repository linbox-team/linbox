/* linbox/blackbox/diagonal-gf2.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 28, 2002.
 *
 * Added parametrization of VectorCategory tags by VectorTraits. See
 * vector-traits.h for more details.
 *
 * ------------------------------------
 *
 * 
 * ========LICENCE========
 * This file is part of the library LinBox.
 * 
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file blackbox/diagonal-gf2.h
 * @ingroup blackbox
 * @brief Random diagonal matrices and diagonal matrices
 * Class especially meant for diagonal precondtionners
 */

#ifndef __LINBOX_diagonal_gf2_H
#define __LINBOX_diagonal_gf2_H

#include "linbox/field/gf2.h"
#include "linbox/blackbox/diagonal.h"

namespace LinBox
{
	template <>
	class Diagonal<GF2, VectorTraits<Vector<GF2>::Dense>::VectorCategory> : public BlackboxArchetype {
	public:

		typedef GF2                       Field;
		typedef Vector<GF2>::Dense        Vector;
		typedef BlackboxArchetype         Blackbox;
		typedef bool                      Element;

		Diagonal (const Field &, const BitVector &y) :
			_v (y)
		{}

		/// The field.
		const Field& field() const
		{
			return *(new GF2());
		}

		Blackbox *clone() const
		{
			return new Diagonal (*this);
		}


		template <class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const
		{
			linbox_check (y.size () == x.size ());
			linbox_check (y.size () == _v.size ());
			typename InVector::const_iterator j1 = x.begin();
			typename OutVector::iterator i = y.begin();
			BitVector::const_iterator j2 = _v.begin();
			for (; i != y.end (); ++i, ++j1, ++j2)
				*i = *j1 & *j2;
			return y;
		}

		Vector& apply (Vector& y, const Vector& x) const
		{
			linbox_check (y.size () == x.size ());
			linbox_check (y.size () == _v.size ());

			BitVector::word_iterator i = y.wordBegin ();
			BitVector::const_word_iterator j1 = x.wordBegin (), j2 = _v.wordBegin ();

			for (; i != y.wordEnd (); ++i, ++j1, ++j2)
				*i = *j1 & *j2;

			return y;
		}

		template <class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const
		{
			return apply (y, x);
		}

		size_t rowdim () const
		{
			return _v.size ();
		}

		size_t coldim () const
		{
			return _v.size ();
		}

		/** Get an entry and store it in the given value
		 * @param x Element in which to store result
		 * @param i Row index
		 * @param j Column index
		 * @return Reference to x
		 */
		Element &getEntry (Element &x, size_t i, size_t j) const
		{
			return (i==j?x=this->_v[i]:x=false);
		}

	private:

		// Bit vector of elements
		BitVector _v;

	}; // template <Field, Vector> class Diagonal<DenseVectorTag>

}

#endif // __LINBOX_diagonal_gf2_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

