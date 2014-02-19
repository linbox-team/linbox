
/* ring/givaro-poly.h
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

#include <iostream>

#include "linbox-config.h"

#include "linbox/integer.h"
#include "givaro/givpoly1.h"

#ifndef LINBOX_ring_givaro_poly_H
#define LINBOX_ring_givaro_poly_H

namespace LinBox
{
	template<class Field>
	class GivaroPolyRandIter;
	
	template<class Domain>
	class GivaroPoly {
		Domain _pd;
	public:
		typedef Domain Domain_t;
		typedef typename Domain::Element Element;
		typedef typename Domain::Type_t Scalar_t;
		
		typedef GivaroPolyRandIter<GivaroPoly<Domain>> RandIter;
		
		Element zero, one, mOne;
		
		GivaroPoly(){}
		
		GivaroPoly(Domain pd) {
			_pd = pd;
			
			_pd.assign(zero, _pd.zero);
			_pd.assign(one, _pd.one);
			_pd.assign(mOne, _pd.mOne);
		}
		
		GivaroPoly(const GivaroPoly &P)
			: _pd(P._pd), zero(P.zero), one(P.one), mOne(P.mOne) {}
		
		GivaroPoly &operator=(const GivaroPoly &F) {
			return *this;
		}
		
		const Domain_t domain() const {
			return _pd;
		}
		
		const typename Domain_t::Domain_t &subdomain() const {
			return _pd.subdomain();
		}
		
		Element &init(Element &x) const {
			_pd.init(x, Givaro::Degree(0), 0);
			return x;
		}
		
		Element &init(Element &x, const Element &y) const {
			_pd.assign(x, y);
			return x;
		}
		
		Element &init(Element &x, const integer &y) const {
			_pd.init(x, Givaro::Degree(0), 0);
			
			integer q = _pd.characteristic();
			
			integer i = 0;
			integer tmp = y;
			
			while (tmp != 0)
			{
				Element xi;
				_pd.init(xi, Givaro::Degree(i), tmp % q);
				_pd.addin(x, xi);
				
				i++;
				tmp /= q;
			}
			
			return x;
		}
		/*
		Element &init(Element &x, const Element &y) const {
			return _pd.assign(x, y);
		}
		*/
		
		integer &convert(integer &x, Element y) const {
			x = 0;
			integer q = _pd.characteristic();
			Givaro::Degree d = _pd.degree(y);
			
			while (d >= 0) {
				integer tmp;
				Scalar_t e;
				
				x *= q;
				_pd.getEntry(e, d, y);
				x += _pd.subdomain().convert(tmp, e);
				
				d--;
			}
			
			return x;
		}
		
		Element &assign(Element &x, const Element y) const {
			return _pd.assign(x, y);
		}
		
		integer &cardinality(integer &c) const {
			return c = -1;
		}
		
		integer &characteristic(integer &c) const {
			return c = _pd.characteristic(c);
		}
		
		bool areEqual(Element x, Element y) const {
			return _pd.areEqual(x, y);
		}
		
		bool isZero(Element x) const {
			return _pd.isZero(x);
		}
		
		bool isOne(Element x) const {
			return _pd.isOne(x);
		}
		
		bool isMOne(Element x) const {
			return _pd.isMOne(x);
		}
		
		std::ostream &write(std::ostream &os) const {
			return os << "polynomials";
		}
		
		std::istream &read(std::istream &is) {
			return is;
		}
		
		std::ostream &write(std::ostream &os, Element x) const {
			return _pd.write(os, x);
		}
		
		std::istream &read(std::istream &is, Element &x) const {
			_pd.read(is, x);
			return is;
		}
		
		Element &add(Element &x, Element y, Element z) const {
			return _pd.add(x, y, z);
		}
		
		Element &sub(Element &x, Element y, Element z) const {
			return _pd.sub(x, y, z);
		}
		
		Element &mul(Element &x, Element y, Element z) const {
			return _pd.mul(x, y, z);
		}
		
		Element &div(Element &x, Element y, Element z) const {
			return _pd.div(x, y, z);
		}
		
		Element &neg(Element &x, Element y) const {
			return _pd.neg(x, y);
		}
		
		Element &inv(Element &x, Element y) const {
			if (_pd.degree(y) == 0 && !_pd.isZero(y)) {
				Scalar_t tmp;
				_pd.subdomain().inv(tmp, _pd.getEntry(tmp, Givaro::Degree(0), y));
				return _pd.assign(x, tmp);
			}
			
			return _pd.assign(x,zero);
		}
		
		Element &axpy(Element &r, Element a, Element x, Element y) const {
			return _pd.axpy(r,a,x,y);
		}
		
		Element &addin(Element &x, Element y) const {
			return _pd.addin(x, y);
		}
		
		Element &subin(Element &x, Element y) const {
			return _pd.subin(x, y);
		}
		
		Element &mulin(Element &x, Element y) const {
			return _pd.mulin(x, y);
		}
		
		Element &divin(Element &x, Element y) const {
			return _pd.divin(x, y);
		}
		
		Element &negin(Element &x) const {
			return _pd.negin(x);
		}
		
		Element &invin(Element &x) const {
			Element tmp;
			inv(tmp, x);
			x = tmp;
			return x;
		}
		
		Element &axpyin(Element &r, Element a, Element x) const {
			return _pd.axpyin(r, a, x);
		}
		
		// PIR Functions
		
		Element &isDivisor(const Element &a, const Element &b) const {
			if (a == zero) return false;
			if (b == zero) return true;
			
			Element tmp;
			mod(tmp, a, b);
			return tmp == zero;
		}
		
		// g = gcd(a,b)
		Element &gcd(Element &g, Element a, Element b) const {
			return _pd.gcd(g,a,b);
		}
		
		// g = gcd(a,b) = a*s + b*t
		Element &xgcd(Element &g, Element &s, Element &t, Element a, Element b) const {
			return _pd.gcd(g,s,t,a,b);
		}
		
		/**
		 * g = gcd(a,b) = a*s + b*t
		 * 1 = s*u + t*v
		 * u = a/g
		 * v = b/g
		 */
		Element &dxgcd(
			Element &g, 
			Element &s,
			Element &t,
			Element &u,
			Element &v,
			Element a,
			Element b) const
		{
			xgcd(g,s,t,a,b);
			
			div(u,a,g);
			div(v,b,g);
			
			return g;
		}
	};
}

#include "randiter/givaro-poly.h"

#endif // LINBOX_ring_givaro_poly_H
