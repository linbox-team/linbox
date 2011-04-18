
/* linbox/randiter/mersenne-twister.C
 * Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura,
 *               1998 Shawn J. Cokus
 *
 * Written by Makoto Matsumoto, Takuji Nishimura, and Shawn J. Cokus
 * Adapted for the LinBox project by Bradford Hovinen
 *
 * ------------------------------------
 *
 * This file is part of LinBox, licensed under the GNU Lesser General Public
 * License. See COPYING for more information.
 *
 * NB: This file is derived from a source file that is licensed under the
 * LGPL. Thus this file *must not* be made into a header file of any kind!
 *
 * ------------------------------------
 *
 * Implementation of the Mersenne twister pseudo-random number generator,
 * crated by Makoto Matsumoto and Takuji Nishimura. Further information on the
 * Mersenne Twister can be found at
 * http://www.math.keio.ac.jp/~matumoto/emt.html
 *
 * This forms the basic underlying algorithm for most psuedo-random number
 * generation in LinBox.
 */

#include "linbox/linbox-config.h"

#include "linbox/randiter/mersenne-twister.h"
#include "linbox/util/debug.h"

namespace LinBox
{

/* transform [0..2^32] -> [0..1] */
static const double doubleTransform = 2.3283064365386962890625e-10;

static const int N = 624;                // length of state vector
static const int M = 397;                // a period parameter
static const uint32 K = 0x9908B0DFU;     // a magic constant

static inline uint32 hiBit (uint32 u)    // mask all but highest   bit of u
	{ return (u) & 0x80000000U; }
static inline uint32 loBit (uint32 u)    // mask all but lowest    bit of u
	{ return (u) & 0x00000001U; }
static inline uint32 loBits (uint32 u)   // mask     the highest   bit of u
	{ return (u) & 0x7FFFFFFFU; }
static inline uint32 mixBits (uint32 u, uint32 v)   // move hi bit of u to hi bit of v
	{ return hiBit (u) | loBits(v); }

MersenneTwister::MersenneTwister (uint32 seed)
	: _state (N + 1), _left (-1)
{
	setSeed (seed);
}

uint32 MersenneTwister::reload ()
{
	register std::vector<uint32>::iterator p0 = _state.begin ();
	register std::vector<uint32>::iterator p2 = _state.begin () + 2;
	register std::vector<uint32>::iterator pM = _state.begin () + M;
	register uint32 s0, s1;
	register int j;

	if (_left < -1)
		setSeed (4357U);

	_left = N - 1, _next = _state.begin () + 1;

	for(s0 = _state[0], s1 = _state[1], j = N - M + 1; --j; s0 = s1, s1 = *p2++)
		*p0++ = *pM++ ^ (mixBits (s0, s1) >> 1) ^ (loBit (s1) ? K : 0U);

	for(pM = _state.begin (), j = M; --j; s0 = s1, s1 = *p2++)
		*p0++ = *pM++ ^ (mixBits (s0, s1) >> 1) ^ (loBit (s1) ? K : 0U);

	s1 = _state[0], *p0 = *pM ^ (mixBits (s0, s1) >> 1) ^ (loBit (s1) ? K : 0U);
	s1 ^= (s1 >> 11);
	s1 ^= (s1 <<  7) & 0x9D2C5680U;
	s1 ^= (s1 << 15) & 0xEFC60000U;
	return(s1 ^ (s1 >> 18));
}


uint32 MersenneTwister::randomInt ()
{
	uint32 y;

	if (--_left < 0)
		return (reload ());

	y  = *_next++;
	y ^= (y >> 11);
	y ^= (y <<  7) & 0x9D2C5680U;
	y ^= (y << 15) & 0xEFC60000U;
	return (y ^ (y >> 18));
}

/* N.B. The following is adapted from Glib 2.2, g_rand_double
 */

double MersenneTwister::randomDouble ()
{
	double retval;

	do {
		/* We set all 52 bits after the point for this, not only the first
		   32. Thats why we need two calls to randomInt */
		retval = randomInt () * doubleTransform;
		retval = (retval + randomInt ()) * doubleTransform;
	} while (retval >= 1.0);

	return retval;
}

/* N.B. The following is adapted from Glib 2.2, g_rand_int_range
 */

uint32 MersenneTwister::randomIntRange (uint32 start, uint32 end) 
{
	linbox_check (end > start);

	uint32 dist = end - start;
	uint32 random;

	/* All tricks doing modulo calculations do not have a perfect
	 * distribution -> We must use the slower way through gdouble for
	 * maximal quality. */
   
	if (dist <= 0x10000L) { /* 2^16 */
		/* This method, which only calls g_rand_int once is only good
		 * for (end - begin) <= 2^16, because we only have 32 bits set
		 * from the one call to g_rand_int (). */

		/* we are using (trans + trans * trans), because g_rand_int only
		 * covers [0..2^32-1] and thus g_rand_int * trans only covers
		 * [0..1-2^-32], but the biggest double < 1 is 1-2^-52. 
		 */

		double double_rand =
			randomInt () * (doubleTransform + doubleTransform * doubleTransform);
      
		random = (uint32) (double_rand * dist);
	} else {
		/* Now we use g_rand_double_range (), which will set 52 bits for
		   us, so that it is safe to round and still get a decent
		   distribution */
		random = (uint32) randomDoubleRange (0, dist);
	}
 
	return start + random;
}

void MersenneTwister::setSeed (uint32 seed) 
{
	//
	// We initialize _state[0..(N-1)] via the generator
	//
	//   x_new = (69069 * x_old) mod 2^32
	//
	// from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
	// _The Art of Computer Programming_, Volume 2, 3rd ed.
	//
	// Notes (SJC): I do not know what the initial state requirements
	// of the Mersenne Twister are, but it seems this seeding generator
	// could be better.  It achieves the maximum period for its modulus
	// (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
	// x_initial can be even, you have sequences like 0, 0, 0, ...;
	// 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
	// 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
	//
	// Even if x_initial is odd, if x_initial is 1 mod 4 then
	//
	//   the          lowest bit of x is always 1,
	//   the  next-to-lowest bit of x is always 0,
	//   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
	//   the 3rd-from-lowest bit of x 4-cycles        ... 0 1 1 0 0 1 1 0 ... ,
	//   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
	//    ...
	//
	// and if x_initial is 3 mod 4 then
	//
	//   the          lowest bit of x is always 1,
	//   the  next-to-lowest bit of x is always 1,
	//   the 2nd-from-lowest bit of x alternates      ... 0 1 0 1 0 1 0 1 ... ,
	//   the 3rd-from-lowest bit of x 4-cycles        ... 0 0 1 1 0 0 1 1 ... ,
	//   the 4th-from-lowest bit of x has the 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
	//    ...
	//
	// The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
	// 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
	// also does well in the dimension 2..5 spectral tests, but it could be
	// better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
	//
	// Note that the random number user does not see the values generated
	// here directly since reloadMT() will always munge them first, so maybe
	// none of all of this matters.  In fact, the seed values made here could
	// even be extra-special desirable if the Mersenne Twister theory says
	// so-- that's why the only change I made is to restrict to odd seeds.
	//

	register uint32 x = (seed | 1U) & 0xFFFFFFFFU;
	register std::vector<uint32>::iterator s = _state.begin ();
	register int j;

	for (_left = 0, *s++ = x, j = N; --j; *s++ = (x *= 69069U) & 0xFFFFFFFFU) ;
}

}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
