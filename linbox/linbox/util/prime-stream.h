/* linbox/randiter/primes.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
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
class PrimeStream
{
    public:

	/** Constructor
	 *
	 * Construct a PrimeStream object.
	 * @param start Starting point; need not be prime itself
	 * @param move_up true if we should move up from the starting point, false otherwise
	 */
	PrimeStream (Integer &start, bool move_up = true)
		: _curr (start), _move_up (move_up) {}

	~PrimeStream () 
	{}

	/** Get the next prime element
	 *
	 * @param a Place to store the next prime element
	 * @return Reference to next prime element
	 */
	Element &next (Element &a) 
	{

		/** LinBox::Integer doesnot support prevprime */
		/*
		if (_move_up == true) {
			nextprime (_curr, _curr);
			a = _curr;
			_curr += 2L;
		} else {
			prevprime (_curr, _curr);
			a = _curr;
			_curr -= 2L;
		}
		*/

		nextprime (_curr, _curr);
		a = _curr;
		_curr += 2L;

		return a;
	}

	/** Operator form for getting the next prime element
	 */
	PrimeStream<Element> &operator >> (Element &a) 
		{ next (a); return *this; }

    private:

	Integer _curr;
	bool    _move_up;
     
}; // class PrimeStream
 
} // namespace LinBox

#endif // __LINBOX_prime_stream_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
