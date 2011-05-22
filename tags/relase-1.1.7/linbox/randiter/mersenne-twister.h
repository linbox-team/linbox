/* linbox/randiter/mersenne-twister.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * This file is part of LinBox, licensed under the GNU Lesser General Public
 * License. See COPYING for more information.
 *
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

#include "linbox/integer.h"

namespace LinBox
{

class MersenneTwister 
{
    public:
	MersenneTwister (uint32 seed = 0);

	uint32 reload ();
	uint32 randomInt ();
	uint32 randomInt () const
	{ return const_cast<MersenneTwister&>(*this).randomInt();}

	uint32 randomIntRange (uint32 start, uint32 end);
	uint32 randomIntRange (uint32 start, uint32 end) const
	{ return const_cast<MersenneTwister&>(*this).randomIntRange(start,end); }

	double randomDouble ();
	double randomDouble ()  const
	{ return const_cast<MersenneTwister&>(*this).randomDouble(); }

	double randomDoubleRange (double start, double end)
		{ return randomDouble () * (end - start) + start; }
	double randomDoubleRange (double start, double end) const
	{ return const_cast<MersenneTwister&>(*this).randomDoubleRange(start,end); }

	void setSeed (uint32 seed);

    private:
	std::vector<uint32>           _state;
	std::vector<uint32>::iterator _next;
	int                           _left;
};
 
}

#ifdef LinBoxSrcOnly
#include <linbox/randiter/mersenne-twister.C>
#endif
#endif // __LINBOX_mersenne_twister_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
