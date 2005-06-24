/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/randiter/generic.h
 * 2004 june, bds and Dan Roche starting from:
 * Copyright (C) 1999-2001 William J Turner,
 *               2002 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
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
 *
 * See COPYING for license information.
 */

#ifndef __GENERIC_RANDITER_H
#define __GENERIC_RANDITER_H

#include <iostream>
#include <vector>

#include "time.h"
#include "linbox/integer.h"
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
	  This is a generator of random field elements that can be used with
	  any field.  It initializes elements using rand().
	  For prime fields with p < 2^32, a near-uniform distrubution can
	  be expected.  For larger fields or non-prime fields, a near-uniform
	  distribution on an unspecified subset of the elements can be expected.
	 */
	template <class Field>
	class GenericRandIter
	{
	    public:

		typedef typename Field::Element Element;

		/** Constructor from field, sampling size, and seed.
		 * The random field element iterator works in the field F, is seeded
		 * by seed, and it returns any one element with probability no more
		 * than 1/min (size, F.characteristic(c)).
		 * A sampling size of zero means to sample from the entire prime subfield.
		 * A seed of zero means to use some arbitrary seed for the generator.
		 * @param F LinBox field in which to do arithmetic
		 * @param size constant integer reference of sample size from which to 
		 *             sample (default = modulus of field)
		 * @param seed constant integer reference from which to seed random number
		 *             generator (default = 0)
		 */
		GenericRandIter (const Field &F, 
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

#ifdef __LINBOX_XMLENABLED
		// XML LinBox::Reader constructor
		GenericRandIter(LinBox::Reader &R) : _F(R.Down(1))
		{
			if(R.haveError()) return;
			R.Up(1);
			if(!R.expectTagName("randiter")) return;
			if(!R.expectAttributeNum("seed", _seed) || !R.expectAttributeNum("size", _size)) return;

			if(_seed == 0) _seed = time(NULL);

			// re-seed the random number generator
			srand(_seed);

			return;

		}
#endif



		GenericRandIter (const GenericRandIter<Field> &R) 
			: _F (R._F), _size (R._size), _seed (R._seed) {}

		~GenericRandIter () {}
    
		GenericRandIter<Field> &operator=(const GenericRandIter<Field> &R)
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
			return (a = ElementEnvelope <Field> (tmp));
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
			std::string s;
			W.setTagName("randiter");
			W.setAttribute("seed", LinBox::Writer::numToString(s, _seed));
			W.setAttribute("size", LinBox::Writer::numToString(s, _size));

			W.addTagChild();
			if(!_F.toTag(W)) return false;
			W.upToParent();

			return true;
		}
#endif


	    private:

		/// Field in which arithmetic is done
		Field _F;

		/// Sampling size
		integer _size;
    
		/// Seed
		long _seed;

	}; // class GenericRandIter
};
#endif //__GENERIC_RANDITER_H
