/* linbox/randiter/mersenne-twister.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
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

 * ------------------------------------
 *
 * Header file for implementation of the Mersenne twister pseudo-random number
 * generator, crated by Makoto Matsumoto and Takuji Nishimura. Further
 * information on the Mersenne Twister can be found at
 * http://www.math.keio.ac.jp/~matsumoto/emt.html [ dead link.  -bds 2007feb]
 *
 * This forms the basic underlying algorithm for most psuedo-random number
 * generation in LinBox.
 *
 * N.B. This module is tested in test-modular under the random number
 * generator tests.
 */

#ifndef __LINBOX_mersenne_twister_H
#define __LINBOX_mersenne_twister_H

#include <vector>

#ifdef __LINBOX_HAVE_STDINT_H
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
// else ??
#endif
#include <stdint.h>
#ifndef INT32_MAX
#error "INT32_MAX is not defined. It should at least be defined in Givaro..."
#endif
#endif


namespace LinBox
{

	class MersenneTwister {
	public:
		MersenneTwister (uint32_t seed = 0);

		uint32_t reload (); // public?
		// random in [0..2^32) ?
		uint32_t randomInt ();
		uint32_t randomInt () const
		{ return const_cast<MersenneTwister&>(*this).randomInt();}

		// random in [start..end).
		uint32_t randomIntRange (uint32_t start, uint32_t end);
		uint32_t randomIntRange (uint32_t start, uint32_t end) const
		{ return const_cast<MersenneTwister&>(*this).randomIntRange(start,end); }

		// random in [0..1] in some sense ?
		double randomDouble ();
		double randomDouble ()  const
		{ return const_cast<MersenneTwister&>(*this).randomDouble(); }

		// random in [start..end] in some sense ?
		double randomDoubleRange (double start, double end)
		{ return randomDouble () * (end - start) + start; }
		double randomDoubleRange (double start, double end) const
		{ return const_cast<MersenneTwister&>(*this).randomDoubleRange(start,end); }

		void setSeed (uint32_t seed);

	private:
		std::vector<uint32_t>           _state;
		std::vector<uint32_t>::iterator _next;
		int                           _left;
	};

}

#if defined(LinBoxSrcOnly) or defined(LinBoxTestOnly)
#include "linbox/randiter/mersenne-twister.C"
#endif
#endif // __LINBOX_mersenne_twister_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

