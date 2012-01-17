/* linbox/randiter/nonzero.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __LINBOX_randiter_nonzero_H
#define __LINBOX_randiter_nonzero_H

#include "linbox/field/archetype.h"
#include "linbox/randiter/archetype.h"
#include "linbox/element/archetype.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/linbox-config.h"

#include <sys/time.h>
#include <stdlib.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <string>

#endif

namespace LinBox
{
	/** Random iterator for nonzero random numbers
	 *
	 * Wraps around an existing random iterator and ensures that the output
	 * is entirely nonzero numbers.
	 **/
	template <class Field, class RandIter = typename Field::RandIter>
	class NonzeroRandIter {
	public:

		typedef typename Field::Element Element;

		NonzeroRandIter (const Field &F, const RandIter &r) :
			_field (F), _r (r)
		{}

		NonzeroRandIter (const NonzeroRandIter& R) :
			_field (R._field), _r (R._r)
		{}

		~NonzeroRandIter()
		{}

		NonzeroRandIter& operator=(const NonzeroRandIter& R)
		{
			if (this != &R) { // guard against self-assignment
				_field = R._field;
				_r = R._r;
			}

			return *this;
		}

		template<class Elt>
		Elt &random (Elt &a)  const
		{
			do _r.random (a); while (_field.isZero (a));
			return a;
		}

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		ElementAbstract &random (ElementAbstract &a)  const
		{
			Element tmp;

			random (tmp);
			return (a = ElementEnvelope <Field> (tmp));
		}

	private:

		Field    _field;
		RandIter _r;

	}; // class NonzeroRandIter

} // namespace LinBox

#endif // __LINBOX_randiter_nonzero_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

