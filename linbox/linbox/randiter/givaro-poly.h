
/* randiter/givaro-poly.h
 * Copyright (C) 2014 Gavin Harrison,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>,
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
 */

#include <givaro/givpoly1.h>
#include <givaro/giv_randiter.h>

#ifndef __LINBOX_randiter_givaro_poly_H
#define __LINBOX_randiter_givaro_poly_H

namespace LinBox
{
	template<class Field>
	class GivaroPolyRandIter {
		typedef typename Field::Domain_t Domain;
		typedef typename Domain::Domain_t SubDomain;
		
		Field _pd;
		Givaro::GIV_randIter<SubDomain,integer> _randIter;
		
		integer _size;
		integer _seed;
	public:
		typedef typename Field::Element Element;
		
		GivaroPolyRandIter(Field pd, 
			const integer& size = 0,
			const integer& seed = 0) 
		{
			_pd = pd;
			_randIter = Givaro::GIV_randIter<SubDomain,integer>(pd.subdomain(), size, seed);
		}
		
		GivaroPolyRandIter(const GivaroPolyRandIter &R)
			: _pd(R._pd), _randIter(R._randIter), _size(R._size), _seed(R._seed) {}
		
		GivaroPolyRandIter &operator=(const GivaroPolyRandIter &R) {
			return *this;
		}
		
		Element &random(Element &a) const {
			return _pd.domain().random(_randIter, a);
		}
		
		Element &random(Element &a, Givaro::Degree d) const {
			return _pd.domain().random(_randIter, a, d);
		}
		
		Element &random(Element &a, long s) const {
			return _pd.domain().random(_randIter, a, s);
		}
	};
}

#endif // __LINBOX_randiter_givaro_poly_H
