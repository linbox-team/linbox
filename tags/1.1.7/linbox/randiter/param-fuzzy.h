/* linbox/randiter/param-fuzzy.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * updated by bds 8/02
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

#ifndef __LINBOX_randiter_param_fuzzy_H
#define __LINBOX_randiter_param_fuzzy_H

#include <iostream>
#include <vector>
#include <time.h>

#include "linbox/integer.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/field/param-fuzzy.h"

namespace LinBox 
{ 

	class ParamFuzzyRandIter
	{
	    public:

		/// element type
		typedef double Element;

		ParamFuzzyRandIter (/*const ParamFuzzy &F, */
				   const integer &size = 0, 
				   const integer &seed = 0)
			: /*_F (F),*/ _size (size), _seed (seed)
		{ 
			/*if (_size == 0) F.cardinality (_size);*/
			if (_seed == 0) _seed = std::time (NULL);    
		}

		ParamFuzzyRandIter (const ParamFuzzy &F,
				    const integer &size = 0, 
				    const integer &seed = 0)
			: _F (F), _size (size), _seed (seed)
		{ 
			if (_size == 0) F.cardinality (_size);
			if (_seed == 0) _seed = std::time (NULL);    
		}

		ParamFuzzyRandIter (const ParamFuzzyRandIter &R) 
			: /*_F (R._F),*/ _size (R._size), _seed (R._seed) {}

		~ParamFuzzyRandIter () {}
    
		ParamFuzzyRandIter &operator=(const ParamFuzzyRandIter &R)
		{
			if (this != &R) { // guard against self-assignment
				_size = R._size;
				_seed = R._seed;
			}

			return *this;
		}
 
		Element &random (Element &a)  const
		{
			// Create new random elements
			if (_size == 0)
				return (a = Element (rand ()));
			else
				return (a = Element (double (rand ())/RAND_MAX)*double (_size));
		}

		ElementAbstract &random (ElementAbstract &a) const
		{
			Element tmp;

			random (tmp);
			return (a = ElementEnvelope <ParamFuzzy> (tmp));
		}

	    private:

		/// Field in which arithmetic is done
		ParamFuzzy _F;

		/// Sampling size
		integer _size;
    
		/// Seed
		integer _seed;

	}; // class ParamFuzzyRandIter : public ParamFuzzyRandIter

} // namespace LinBox 

#endif // __LINBOX_randiter_param_fuzzy_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
