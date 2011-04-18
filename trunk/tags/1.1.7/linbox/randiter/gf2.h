/* linbox/randiter/gf2.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_randiter_gf2_H
#define __LINBOX_randiter_gf2_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
#include "linbox/field/modular.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/util/commentator.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <string>

#endif

namespace LinBox 
{ 

class GF2RandIter
{
    public:
	typedef bool Element;

	/** Constructor from field, sampling size, and seed.
	 * The random field element iterator works in the field F, is seeded
	 * by seed, and it returns any one element with probability no more
	 * than 1/min (size, F.cardinality (c)).
	 * A sampling size of zero means to sample from the entire field.
	 * A seed of zero means to use some arbitrary seed for the generator.
	 * Purely virtual.
	 * @param F LinBox field archetype object in which to do arithmetic
	 * @param size constant integer reference of sample size from which to 
	 *             sample (default = modulus of field)
	 * @param seed constant integer reference from which to seed random number
	 *             generator (default = 0)
	 */
	GF2RandIter (const GF2 &, 
		     const integer &  = 0 , 
		     const integer &seed = 0)
	{
		long _seed = seed;

		if (_seed == 0) _seed = time (NULL);
		MT.setSeed (_seed);
	}

	GF2RandIter (const GF2RandIter &) {}

	/** Destructor.
	 * This destructs the random field element generator object.
	 */
	~GF2RandIter () {}
    
	/** Assignment operator.
	 * Assigns ModularRandIter object R to generator.
	 * @param  R ModularRandIter object.
	 */
	GF2RandIter &operator = (const GF2RandIter &)
		{ return *this; }
 
	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	bool &random (bool &a)  const
		{ return a = MT.randomIntRange (0, 2); }

	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	BitVector::reference random (BitVector::reference a)  const
		{ return a = MT.randomIntRange (0, 2); }


	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	std::_Bit_reference random (std::_Bit_reference a)  const
		{ return a = MT.randomIntRange (0, 2); }

	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	ElementAbstract &random (ElementAbstract &a)  const
	{
		bool tmp;

		random (tmp);
		return (a = ElementEnvelope <GF2> (tmp));
	}


	uint32& random (uint32& a)  const
		{ return a = MT.randomInt(); }

    MersenneTwister& getMT() { return MT; }
    const MersenneTwister& getMT() const { return MT; }

    private:

	MersenneTwister MT;

}; // class GF2RandIter

} // namespace LinBox 

#endif // __LINBOX_randiter_gf2_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
