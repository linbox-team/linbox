/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/randiter/nonzero.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __RANDITER_NONZERO_H
#define __RANDITER_NONZERO_H

#include "linbox/field/archetype.h"
#include "linbox/randiter/archetype.h"
#include "linbox/element/archetype.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox-config.h"

#include <sys/time.h>
#include <stdlib.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

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
	class NonzeroRandIter
	{
	    public:
    
		typedef typename Field::Element Element;

		NonzeroRandIter (const Field &F, const RandIter &r)
			: _F (F), _r (r) {}

		NonzeroRandIter (const NonzeroRandIter& R)
			: _F (R._F), _r (R._r) {}

#ifdef __LINBOX_XMLENABLED
		NonzeroRandIter(Reader &R) : _F(R.Down(1)), _r(R.Up(1)) {}
#endif

		~NonzeroRandIter() 
			{}

		NonzeroRandIter& operator=(const NonzeroRandIter& R)
		{
			if (this != &R) { // guard against self-assignment
				_F = R._F;
				_r = R._r;
			}

			return *this;
		}

		Element &random (Element &a)  const
		{
			do _r.random (a); while (_F.isZero (a));
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

#ifdef __LINBOX_XMLENABLED
		ostream &write(ostream &os) const 
		{
			return _r.write(os);
		}

		bool toTag(Writer &W) const
		{
			return _r.toTag(W);
		}
#endif


	    private:

		Field    _F;
		RandIter _r;
     
	}; // class NonzeroRandIter
 
} // namespace LinBox

#endif // __RANDITER_NONZERO_H