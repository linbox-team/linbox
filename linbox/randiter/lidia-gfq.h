/* Copyright (C) 2010 LinBox
 * Written by  Pascal Giorgi
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/* File: src/wrapper/by_scope/field/LIDIA_randiter.h
 * Author: Pascal Giorgi for the LinBox group
 */

#ifndef __LINBOX_lidia_randiter_gfq_H
#define __LINBOX_lidia_randiter_gfq_H

#include "LiDIA/gf_element.h"

#include "linbox/integer.h"
#include "linbox/linbox-config.h"

namespace LinBox
{

	template<class field> class LidiaGfqRandIter
	{
	public:


		/* Element of the class Field
		*/
		typedef typename field::Element Element;

		/* Constructor of the class LidiaGfqRandIter .
		 * the size has to be a divisor of the degree of the Field F.
		 * the seed doesn't do something for the moment.
		 */
		LidiaGfqRandIter(const field& F,
				 const integer& size = 0,
				 const uint64_t seed = 0) :
			_size(size), _seed(seed) , GF(F)
		{
                    if (_seed == 0) _seed = static_cast<uint64_t>(std::time(nullptr));

			integer cardinality ;
			F.cardinality(cardinality);
			if ( (cardinality != integer(-1)) &&  (_size > cardinality) )
				_size = cardinality;


#ifdef TRACE
			cout << "created random generator with size " << _size
			<< " and seed " << _seed << endl;
#endif // TRACE

			// Seed random number generator
			//srand(_seed);
			LiDIA::bigint tmp;

			string_to_bigint((std::string(seed)).c_str(),tmp);
			LiDIA::bigint::seed(tmp);


		}

		LidiaGfqRandIter(const LidiaGfqRandIter& R) :
			_size(R._size), _seed(R._seed), GF(R.GF) {}


		~LidiaGfqRandIter(void) {}


		LidiaGfqRandIter& operator=(const LidiaGfqRandIter& R)
		{
			if (this != &R) // guard against self-assignment
			{
				_size = R._size;
				_seed = R._seed;
			}

			return *this;
		}


		Element& random (Element& x)  const
		{
			Element e(GF);
			if (_size == 0)
				e.randomize();
			else
				e.randomize(_size);

			x.assign(e);
			return x;

			/* Another method allowing generate with a sampling size

			   long temp;
			   Element e;
			   if (_size==0)
			   temp= rand();
			   else
			   temp= static_cast<long>((double(rand())/RAND_MAX)*double(_size));

			   GF.init(e,integer(temp));
			   return x=e;
			   */

		}

		LidiaGfqRandIter(void) :
		       	_size(0), _seed(static_cast<uint64_t>(std::time(nullptr)))
		{}

	private:

		/// Sampling size
		integer _size;

		/// Seed
		uint64_t _seed;

		field GF;

	}; // class LidiaGfqRandIter

} // namespace LinBox

#endif // __LINBOX_lidia_randiter_gfq_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
