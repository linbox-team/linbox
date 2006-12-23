/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/randiter/modular.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Changed LargeModularRandIter to ModularRandIter, parameterized on the
 * element type. This change is for compatibility with the changes in
 * field/modular.h
 *
 * Renamed from large-modular.h to modular.h
 * ------------------------------------
 * 2002-05-14 William J. Turner <wjturner@acm.org>
 *
 * Seeded random number generator in constructor.  _seed was never used
 * before.
 * ------------------------------------
 * 2005-06-24 William J. Turner <wjturner@acm.org>
 *
 * Removed using declarations.
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __RANDITER_MODULAR_H
#define __RANDITER_MODULAR_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
#include "linbox/field/modular.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/util/commentator.h"
#include "linbox/randiter/mersenne-twister.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

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
	template <class Element>
	class ModularRandIter
	{
	    public:

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
		ModularRandIter (const Modular<Element> &F, 
				 const integer &size = 0, 
				 const integer &seed = 0)
			: _F (F), _size (size), _seed (seed)
		{
			if (_seed == 0) _seed = time (NULL);

			integer cardinality;

			F.cardinality (cardinality);

			if ((_size == 0) || (_size > cardinality))
				_size = cardinality;

			commentator.report (10, INTERNAL_DESCRIPTION)
				<< "Created random generator with size " << _size 
				<< " and seed " << _seed << std::endl;

			// Seed random number generator
			srand (_seed);
		}

		/** Copy constructor.
		 * Constructs ModularRandIter object by copying the random field
		 * element generator.
		 * This is required to allow generator objects to be passed by value
		 * into functions.
		 * @param  R ModularRandIter object.
		 */
		ModularRandIter (const ModularRandIter<Element> &R) 
			: _F (R._F), _size (R._size), _seed (R._seed) {}

		/** Destructor.
		 * This destructs the random field element generator object.
		 */
		~ModularRandIter () {}
    
		/** Assignment operator.
		 * Assigns ModularRandIter object R to generator.
		 * @param  R ModularRandIter object.
		 */
		ModularRandIter<Element> &operator=(const ModularRandIter<Element> &R)
		{
			if (this != &R) { // guard against self-assignment
				_size = R._size;
				_seed = R._seed;
			}

			return *this;
		}
 
		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		Element &random (Element &a) const
		{ return _F.init(a,rand()); }

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * Required by abstract base class.
		 * @return reference to random field element
		 */
		ElementAbstract &random (ElementAbstract &a) const
		{
			Element tmp;

			random (tmp);
			return (a = ElementEnvelope <Modular<Element> > (tmp));
		}

	    private:

		/// Field in which arithmetic is done
		Modular<Element> _F;

		/// Sampling size
		integer _size;
    
		/// Seed
		long _seed;

	}; // class ModularRandIter

	template <class Element>
	class ModularBase<Element>::RandIter {
		ModularRandIter<Element> _r;

	    public:
		RandIter (const Modular<Element> &F, const integer &size = 0, const integer &seed = 0)
			: _r (F, size, seed) {}
		RandIter (const ModularBase<Element>::RandIter &r)
			: _r (r._r) {}

		~RandIter () {}
		RandIter &operator= (const RandIter &r)
			{ _r = r._r; return *this; }
		Element &random (Element &a) const
			{ return _r.random (a); }
		ElementAbstract &random (ElementAbstract &a) const
			{ return _r.random (a); }

	};

	template <>
	class ModularBase<uint16>::RandIter {
		MersenneTwister _r;
		uint16 _size;
		uint16 _seed;

	    public:
		typedef uint16 Element;

		RandIter (const Modular<Element> &F, const integer &size = 0, const integer &seed = 0)
		{
			_seed = seed;
			_size = size;

			if (_seed == 0) _seed = time (NULL);

			integer c;

			F.cardinality (c);

			linbox_check (c != -1);

			if ((_size == 0) || (_size > double (c)))
				_size = c;

			_r.setSeed (_seed);
		}

		RandIter (const ModularBase<Element>::RandIter &r)
			: _r (r._r), _size (r._size), _seed (r._seed) {}
#ifdef __LINBOX_XMLENABLED
		RandIter(LinBox::Reader &R)
		{
			if(!R.expectTagName("randiter")) return;
			if(!R.expectAttributeNum("seed",_seed) || !R.expectAttributeNum("size", _size)) return;

			if(_seed == 0) _seed = time(NULL);

			_r.setSeed(_seed);
			
			return;
		}
#endif


		~RandIter () {}
		RandIter &operator= (const RandIter &r)
			{ _r = r._r; return *this; }
		Element &random (Element &a) const
			{ return a = _r.randomIntRange (0, _size); }
		ElementAbstract &random (ElementAbstract &a)  const
			{ return a = ElementEnvelope <Modular<Element> >
				  (_r.randomIntRange (0, _size)); }

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
			std::string s;
			W.setTagName("randiter");
			W.setAttribute("seed", LinBox::Writer::numToString(s, _seed));
			W.setAttribute("size", LinBox::Writer::numToString(s, _size));

			return true;
		}
#endif		


	};

	template <>
	class ModularBase<uint32>::RandIter {
		MersenneTwister _r;
		uint32 _size;
		uint32 _seed;

	    public:
		typedef uint32 Element;

		RandIter (const Modular<Element> &F, const integer &size = 0, const integer &seed = 0)
		{
			_seed = seed;
			_size = size;

			if (_seed == 0) _seed = time (NULL);

			integer c;

			F.cardinality (c);

			linbox_check (c != -1);

			if ((_size == 0) || (_size > double (c)))
				_size = c;

			_r.setSeed (_seed);
		}

		RandIter (const ModularBase<Element>::RandIter &r)
			: _r (r._r), _size (r._size), _seed (r._seed) {}
#ifdef __LINBOX_XMLENABLED
		RandIter(LinBox::Reader &R) {
			if(!R.expectTagName("randiter")) return;
			if(!R.expectAttributeNum("seed", _seed) || !R.expectAttributeNum("size", _size)) return;

			if(_seed == 0) _seed = time(NULL);

			_r.setSeed(_seed);

			return;
		}
#endif

		~RandIter () {}
		RandIter &operator= (const RandIter &r)
			{ _r = r._r; return *this; }
		Element &random (Element &a) const
			{ return a = _r.randomIntRange (0, _size); }
		ElementAbstract &random (ElementAbstract &a) const
			{ return a = ElementEnvelope <Modular<Element> >
				  (_r.randomIntRange (0, _size)); }
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
			std::string s;

			W.setTagName("randiter");
			W.setAttribute("seed", LinBox::Writer::numToString(s, _seed));
			W.setAttribute("size", LinBox::Writer::numToString(s, _size));

			return true;
		}
#endif
				
			

	};
} // namespace LinBox 

#endif // _LARGE_MODULAR_RANDITER_
