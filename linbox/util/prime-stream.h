/* linbox/randiter/primes.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
 *
 * Stream of prime numbers
 */

#ifndef __LINBOX_prime_stream_H
#define __LINBOX_prime_stream_H

#include "linbox/integer.h"

#include <sys/time.h>
#include <stdlib.h>

namespace LinBox
{

	/** Prime number stream
	 *
	 * Provides a source of prime numbers of given characteristics to use in parts
	 * of the library that need to reduce modulo one or more primes.
	 **/
	template <class Element>
	class PrimeStream {
	public:

		/** Constructor
		 *
		 * Construct a PrimeStream object.
		 * @param start Starting point; need not be prime itself
		 * @param move_up true if we should move up from the starting point, false otherwise
		 */
		PrimeStream (Integer &start, bool move_up = true) :
			_curr (start), _move_up (move_up)
		{}

		~PrimeStream ()
		{}

		/** Get the next prime element
		 *
		 * @param a Place to store the next prime element
		 * @return Reference to next prime element
		 */
		Element &next (Element &a)
		{

			if (_move_up == true) {
				_IPD.nextprimein(_curr);
				a = _curr;
				_curr += 2L;
			}
			else {
				_IPD.prevprimein(_curr);
				a = _curr;
				_curr -= 2L;
			}

			return a;
		}

		/** Operator form for getting the next prime element
		*/
		PrimeStream<Element> &operator >> (Element &a)
		{ next (a); return *this; }

	private:

		Integer _curr;
		bool    _move_up;
		Givaro::IntPrimeDom _IPD; //!< empty struct dealing with primality

	}; // class PrimeStream

} // namespace LinBox

#endif // __LINBOX_prime_stream_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
