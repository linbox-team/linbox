/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/randiter/gmp-rational.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __RANDITER_GMP_RATIONAL_H
#define __RANDITER_GMP_RATIONAL_H

#include "linbox/field/gmp-rational.h"
#include "linbox/element/gmp-rational.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>

using std::ostream;
using std::string;

#endif

#include <sys/time.h>
#include <stdlib.h>

namespace LinBox
{

class GMPRationalRandIter
{
    public:
    
	typedef GMPRationalElement Element;
    
	GMPRationalRandIter (const GMPRationalField &F,
			     const integer &size = 0,
			     const integer &seed = 0)
		: _F (F), _size (size), _seed (seed)
	{
		if (seed == 0)
			_seed = time (NULL);
	}



	GMPRationalRandIter (const GMPRationalRandIter& R)
		: _F (R._F), _size (R._size), _seed (R._seed) {}


#ifdef __LINBOX_XMLENABLED
	GMPRationalRandIter(Reader &R) : _F(R.Down(1))
	{
		R.Up(1);
		if(!R.expectTagName("randiter")) return;
		if(!R.expectAttributeNum("seed", _seed) || !R.expectAttributeNum("size", _size)) return;

		if(_seed == 0) _seed = time( NULL);

		return;

	}
#endif


	~GMPRationalRandIter() 
	{}
    
	GMPRationalRandIter& operator=(const GMPRationalRandIter& R)
	{
		if (this != &R) { // guard against self-assignment
			_F = R._F;
			_seed = R._seed;
			_size = R._size;
		}
		return *this;
	}
 
	Element &random (Element &a)  const
	{
		unsigned int s;
		int value = 0;

		if (_size == 0) {
			s = _seed;

			value = rand_r (&s);

			mpz_set_si (mpq_numref (a.rep), value);

			do {
				value = rand_r (&s);
			} while (value == 0);

			const_cast<integer&>(_seed) = s;
			mpz_set_si (mpq_denref (a.rep), value);
		}
		else {
			unsigned int s;
			int num, den;

			s = _seed;
			num = rand_r (&s);

			if (_size > 0) {
				unsigned long tmp = _size;
				num %= tmp;
				den = 1L;
			} else {
				den = rand_r (&s);
			}

			const_cast<integer&>(_seed) = s;

			mpz_set_si (mpq_numref (a.rep), num);
			mpz_set_si (mpq_denref (a.rep), den);
		}

		mpq_canonicalize (a.rep);

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
		return (a = ElementEnvelope <GMPRationalField> (tmp));
	}

#ifdef __LINBOX_XMLENABLED
	ostream &write(ostream &os) const
	{
		Writer W;
		if( toTag(W))
			W.write(os);

		return os;
	}

	bool toTag(Writer &W) const
	{
		string s;
		W.setTagName("randiter");
		W.setAttribute("seed", Writer::numToString(s, _seed));
		W.setAttribute("size", Writer::numToString(s, _size));

		W.addTagChild();
		if(!_F.toTag(W)) return false;
		W.upToParent();

		return true;
	}
#endif

    private:

	GMPRationalField _F;

	integer _size;
	integer _seed;
     
}; // class GMPRationalRandIter
 
} // namespace LinBox

#endif // __RANDITER_GMP_RANDOM_H
