
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

#include "linbox/integer.h"
#include <givaro/givpoly1.h>
#include <givaro/givquotientdomain.h>

#ifndef LINBOX_ring_givaro_poly_quotient_H
#define LINBOX_ring_givaro_poly_quotient_H

namespace LinBox
{
	template<class Field>
	class GivaroPolyRandIter;

	template<class Domain>
	class GivaroPolyQuotient : public Givaro::QuotientDom<Domain> {
		Domain _pd;
	public:
		typedef Domain Domain_t;
		typedef typename Domain::Element Element;
		typedef typename Domain::Type_t Scalar_t;

		typedef GivaroPolyRandIter<GivaroPoly<Domain> > RandIter;

		Element zero, one, mOne, _f;

		GivaroPolyQuotient(){}

		GivaroPolyQuotient(const Domain &pd, Element f) : 
			Givaro::QuotientDom<Domain>(pd, f), _pd(pd) {
			_pd.assign(zero, _pd.zero);
			_pd.assign(one, _pd.one);
			_pd.assign(mOne, _pd.mOne);
			_pd.assign(_f, f);
		}

		GivaroPolyQuotient(const GivaroPolyQuotient &P)
			: Givaro::QuotientDom<Domain>(P._pd, P._f), _pd(P._pd), 
			zero(P.zero), one(P.one), mOne(P.mOne), _f(_pd._f) {}

		GivaroPolyQuotient &operator=(const GivaroPolyQuotient &F) {
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

		Element &init(Element &x, const integer &y) const {
			_pd.init(x, Givaro::Degree(0), 0);

			integer q = subdomain().cardinality();

			integer i = 0;
			integer tmp = y;

			while (tmp != 0)
			{
				Element xi;
				_pd.init(xi, Givaro::Degree(int64_t(i)), tmp % q);
				_pd.addin(x, xi);

				i++;
				tmp /= q;
			}

			return _pd.modin(x, _f);
		}
		
		Element &init(Element &x, const std::vector<integer> &y) const {
			_pd.init(x, Givaro::Degree(0), 0);
			
			for (int i = 0; i < y.size(); i++)
			{
				Element xi;
				_pd.init(xi, Givaro::Degree(i), y[i]);
				_pd.addin(x, xi);
			}
			
			return _pd.modin(x, _f);
		}

		integer &convert(integer &x, const Element &y) const {
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

		Element &assign(Element &x, const Element &y) const {
			return _pd.assign(x, y);
		}

		integer &cardinality(integer &c) const {
			integer tmp;
			return c = _pd.degree(_f).value() * _pd.characteristic(tmp);
		}

		integer &characteristic(integer &c) const {
			return _pd.characteristic(c);
		}

		//bool isZero(const Element &x) const {
		//	return _pd.isZero(x);
		//}

		//bool isOne(const Element &x) const {
		//	return _pd.isOne(x);
		//}

		//bool isMOne(const Element &x) const {
		//	return _pd.isMOne(x);
		//}

		//bool areEqual(const Element &x, const Element &y) const {
		//	return _pd.areEqual(x, y);
		//}
		
		Element &normalize(Element &z, const Element &x) const {
			return _pd.gcd(z, x, _f);
		}

		Element &normalizeIn(Element &x) const {
			return _pd.gcdin(x, _f);
		}
		
		bool areAssociates(const Element &x, const Element &y) const {
			Element a, b;
			_pd.gcd(a, x, _f);
			_pd.gcd(b, y, _f);
			return _pd.degree(a).value() == _pd.degree(b).value();
		}

		//std::ostream &write(std::ostream &os) const {
		//	return os << "polynomials";
		//}

		//std::istream &read(std::istream &is) {
		//	return is;
		//}

		//std::ostream &write(std::ostream &os, const Element &x) const {
		//	return _pd.write(os, x);
		//}

		//std::istream &read(std::istream &is, Element &x) const {
		//	_pd.read(is, x);
		//	return is;
		//}

		Element &add(Element &x, const Element &y, const Element &z) const {
			return _pd.add(x, y, z);
		}

		Element &sub(Element &x, const Element &y, const Element &z) const {
			return _pd.sub(x, y, z);
		}

		//Element &mul(Element &x, const Element &y, const Element &z) const {
		//	return _pd.mul(x, y, z);
		//}

		//Element &div(Element &x, const Element &y, const Element &z) const {
		//	return _pd.div(x, y, z);
		//}

		Element &neg(Element &x, const Element &y) const {
			return _pd.neg(x, y);
		}

		//Element &inv(Element &x, const Element &y) const {
		//	if (_pd.degree(y) == 0 && !_pd.isZero(y)) {
		//		Scalar_t tmp;
		//		_pd.subdomain().inv(tmp, _pd.getEntry(tmp, Givaro::Degree(0), y));
		//		return _pd.assign(x, tmp);
		//	}
		//	return _pd.assign(x,zero);
		//}

		//Element &axpy(Element &r, const Element &a, const Element &x, const Element &y) const {
		//	return _pd.axpy(r,a,x,y);
		//}

		//Element &addin(Element &x, const Element &y) const {
		//	return _pd.addin(x, y);
		//}

		//Element &subin(Element &x, const Element &y) const {
		//	return _pd.subin(x, y);
		//}

		//Element &mulin(Element &x, const Element &y) const {
		//	return _pd.mulin(x, y);
		//}

		//Element &divin(Element &x, const Element &y) const {
		//	return _pd.divin(x, y);
		//}

		//Element &negin(Element &x) const {
		//	return _pd.negin(x);
		//}

		//Element &invin(Element &x) const {
		//	Element tmp;
		//	inv(tmp, x);
		//	x = tmp;
		//	return x;
		//}

		//Element &axpyin(Element &r, const Element &a, const Element &x) const {
		//	return _pd.axpyin(r, a, x);
		//}

		// PIR Functions

		bool isDivisor(const Element &x, const Element &y) const {
			Element a, b;
			_pd.gcd(a, x, _f);
			_pd.gcd(b, y, _f);
			return _pd.degree(a).value() >= _pd.degree(b).value();
		}
		
		// a = q b + r
		//Element &quo(Element &q, const Element &a, const Element &b) const {
		//	return div(q,a,b);
		//}
		
		//Element &rem(Element &r, const Element &a, const Element &b) const {
		//	return _pd.mod(r,a,b);
		//}
		
		//Element &divrem(Element &q, Element &r, const Element &a, const Element &b) const {
		//	return _pd.divmod(q,r,a,b);
		//}
		
		//Element &modin(Element &a, const Element &b) const {
		//	return _pd.modin(a, b);
		//}

		//Element& normalIn(Element &x) const {
		//	return _pd.gcdin(x, _f);
		//}
		
		//bool isUnit(const Element& x) const {
		//	Element tmp;
		//	return _pd.degree(_pd.gcd(tmp, x, _f)).value() == 0;
		//}

		// g = gcd(a,b)
		Element& gcd(Element &g, const Element &a, const Element &b) const {
			Element ga, gb;
			_pd.gcd(ga,a,_f);
			_pd.gcd(gb,b,_f);
			
			if (_pd.degree(ga).value() <= _pd.degree(gb).value()) {
				return _pd.assign(g, ga);
			} else {
				return _pd.assign(g, gb);
			}
		}
		
		Element& gcdin(Element &a, const Element &b) const {
			Element ga, gb;
			_pd.gcd(ga,a,_f);
			_pd.gcd(gb,b,_f);
			
			if (_pd.degree(ga).value() <= _pd.degree(gb).value()) {
				return _pd.assign(a, ga);
			} else {
				return _pd.assign(a, gb);
			}
		}
		
		// g = gcd(a,b)
		//Element& gcd(Element &g, Element &s, Element &t, const Element &a, const Element &b) const {
		//	return _pd.gcd(g,s,t,a,b);
		//}

		// g = gcd(a,b) = a*s + b*t
		//Element &xgcd(Element &g, Element &s, Element &t, const Element &a, const Element &b) const {
		//	return _pd.gcd(g,s,t,a,b);
		//}

		/**
		 * g = gcd(a,b) = a*s + b*t
		 * 1 = s*u + t*v
		 * u = a/g
		 * v = b/g
		 */
		//Element &dxgcd(
		//	Element &g,
		//	Element &s,
		//	Element &t,
		//	Element &u,
		//	Element &v,
		//	const Element &a,
		//	const Element &b) const
		//{
		//	xgcd(g,s,t,a,b);

		//	div(u,a,g);
		//	div(v,b,g);

		//	return g;
		//}
	};
}

#include "linbox/randiter/givaro-poly.h"

#endif // LINBOX_ring_givaro_poly_quotient_H
