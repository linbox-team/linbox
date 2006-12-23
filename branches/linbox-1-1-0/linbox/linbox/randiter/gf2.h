/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/randiter/gf2.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __RANDITER_GF2_H
#define __RANDITER_GF2_H

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
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <string>

#endif

namespace LinBox 
{ 

/** Random field base element generator.
 * This encapsulated class is a generator of random field base elements for 
 * the encapsulating field.
 * It is required to contain constructors from a field object and
 * two integers.  The first integer being a cardinality of a set to 
 * draw the random elements from, and the second being a seed for the 
 * random number generator.
 * It is also required to contain a copy constructor, a destructor, and
 * an operator () which acts on a reference to a field base element.  In this 
 * operator (), the random element is placed into the input field base element 
 * and also returned as a reference.
 */

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
	GF2RandIter (const GF2 &F, 
		     const integer &size = 0, 
		     const integer &seed = 0)
	{
		long _seed = seed;

		if (_seed == 0) _seed = time (NULL);
		_MT.setSeed (_seed);
	}

	/** Copy constructor.
	 * Constructs ModularRandIter object by copying the random field
	 * element generator.
	 * This is required to allow generator objects to be passed by value
	 * into functions.
	 * @param  R ModularRandIter object.
	 */
	GF2RandIter (const GF2RandIter &R) {}

#ifdef __LINBOX_XMLENABLED
	GF2RandIter(LinBox::Reader &R)
	{
		long seed, size;
		R.Up(1);

		if(!R.expectTagName("randiter")) return;
		if(!R.expectAttributeNum("seed", seed) || !R.expectAttributeNum("size", size)) return;

		if(size != 2) {
			R.setErrorString("Got GF(2) randiter w/ size not 2");
			R.setErrorCode(LinBox::Reader::OTHER);
			return;
		}
		if(seed == 0) seed == time(NULL);
		_MT.setSeed(seed);

		return;
	}
#endif	       

	


	/** Destructor.
	 * This destructs the random field element generator object.
	 */
	~GF2RandIter () {}
    
	/** Assignment operator.
	 * Assigns ModularRandIter object R to generator.
	 * @param  R ModularRandIter object.
	 */
	GF2RandIter &operator = (const GF2RandIter &R)
		{ return *this; }
 
	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	bool &random (bool &a)  const
		{ return a = _MT.randomIntRange (0, 2); }

	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	BitVector::reference random (BitVector::reference a)  const
		{ return a = _MT.randomIntRange (0, 2); }


	/** Random field element creator.
	 * This returns a random field element from the information supplied
	 * at the creation of the generator.
	 * Required by abstract base class.
	 * @return reference to random field element
	 */
	std::_Bit_reference random (std::_Bit_reference a)  const
		{ return a = _MT.randomIntRange (0, 2); }

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

#ifdef __LINBOX_XMLENABLED
	std::ostream &write(std::ostream &os) const
	{
		LinBox::Writer W;
		if( toTag(W))
			W.write(os);

		return os;
	}

	bool toTag(LinBox::Writer &W) const
	{
		W.setTagName("randiter");
		W.setAttribute("seed", "0");
		W.setAttribute("size", "2");
		
		return true;
	}
#endif


    private:

	MersenneTwister _MT;

}; // class GF2RandIter

} // namespace LinBox 

#endif // __RANDITER_GF2_H
