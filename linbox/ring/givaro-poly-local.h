
/* ring/givaro-poly-local.h
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

#include <linbox/integer.h>
#include <linbox/polynomial/dense-polynomial.h>
#include <linbox/ring/polynomial-ring.h>

#ifndef LINBOX_ring_givaro_poly_local_H
#define LINBOX_ring_givaro_poly_local_H

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
                                   const integer& seed = 0) :
                _randIter(Givaro::GIV_randIter<SubDomain,integer>(pd.subdomain(), uint64_t(size), seed))
                {_pd = pd;}

                GivaroPolyRandIter(const GivaroPolyRandIter &R)
                        : _pd(R._pd), _randIter(R._randIter), _size(R._size), _seed(R._seed) {}

                GivaroPolyRandIter &operator=(const GivaroPolyRandIter &R) {
                        return *this;
                }

                Element &random(Element &a)  {
                        return _pd.domain().random(_randIter, a);
                }

                Element &random(Element &a, Givaro::Degree d)  {
                        return _pd.domain().random(_randIter, a, d);
                }

                Element &random(Element &a, long s) const {
                        return _pd.domain().random(_randIter, a, s);
                }
        };
}


namespace LinBox
{
	template<class Field>
	class GivaroPolyRandIter;

	template<class Domain>
	class GivaroPolyLocal {
		Domain _pd;
	public:
		typedef Domain Domain_t;
		typedef typename Domain_t::Type_t Scalar_t;
		typedef typename Domain_t::Element Polynomial_t;
		
		typedef struct {
			Polynomial_t f;
			size_t e;
		} Element;

		typedef GivaroPolyRandIter<PolynomialRing<Domain> > RandIter;

		Element zero, one, mOne;
		Polynomial_t _f;
		size_t _e;
		Polynomial_t _fe;

		GivaroPolyLocal(){}

		GivaroPolyLocal(
			const Domain &pd,
			const Polynomial_t f,
			const size_t e
		) : _pd(pd), _f(f), _e(e)
		{
			_pd.pow(_fe, f, e);
			zero = { _pd.zero, 0 };
			one = { _pd.one, 0 };
			mOne = { _pd.mOne, 0 };
		}

		GivaroPolyLocal(const GivaroPolyLocal &P) : _pd(P._pd), _e(P._e)
		{
			_pd.assign(_f, P._f);
			_pd.assign(_fe, P._fe);
			assign(zero, P.zero);
			assign(one, P.one);
			assign(mOne, P.mOne);
		}

		GivaroPolyLocal &operator=(const GivaroPolyLocal &F) {
			return *this;
		}

		const Domain_t domain() const {
			return _pd;
		}

		const typename Domain_t::Domain_t &subdomain() const {
			return _pd.subdomain();
		}

		Element &init(Element &x) const {
			_pd.assign(x.f, zero);
			x.e = 0;
			return x;
		}
		
		Element &init(Element &x, const Polynomial_t y) const {			
			Polynomial_t g, a;
			_pd.mod(a, y, _fe);
			_pd.gcd(g, a, _fe);
			x.e = _pd.degree(g).value() / _pd.degree(_f).value();
			
			if (x.e >= _e) {
				return assign(x, zero);
			}
			
			_pd.div(x.f, a, _pd.pow(g, _f, x.e));
			return x;
		}

		Element &init(Element &x, const integer &y) const {
			Polynomial_t a;
			_pd.init(a, Givaro::Degree(0), 0);

			integer q = subdomain().cardinality();

			integer i = 0;
			integer tmp = y;

			while (tmp != 0)
			{
				Polynomial_t xi;
				_pd.init(xi, Givaro::Degree(int64_t(i)), tmp % q);
				_pd.addin(a, xi);

				i++;
				tmp /= q;
			}

			return init(x, a);
		}
		
		Element &init(Element &x, const std::vector<integer> &y) const {
			Polynomial_t a;
			
			_pd.init(a, Givaro::Degree(0), 0);
			
			for (size_t i = 0; i < y.size(); i++)
			{
				Element xi;
				_pd.init(xi, Givaro::Degree(i), y[i]);
				_pd.addin(a, xi);
			}
			
			return init(x, a);
		}
		
		Polynomial_t &convert(Polynomial_t &x, const Element &y) const {
			_pd.pow(x, _f, y.e);
			_pd.mulin(x, y.f);
			return _pd.modin(x, _fe);
		}

		integer &convert(integer &x, const Element &z) const {
			Polynomial_t y;
			_pd.convert(y, z);
			
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
			_pd.assign(x.f, y.f);
			x.e = y.e;
			return x;
		}

		integer &cardinality(integer &c) const {
			integer tmp;
			return c = _pd.degree(_fe).value() * _pd.characteristic(tmp);
		}

		integer &characteristic(integer &c) const {
			return c = _pd.characteristic(c);
		}

		bool areEqual(const Element &x, const Element &y) const {
			return x.e == y.e && _pd.areEqual(x.f, y.f);
		}
		
		bool areEqual(const Element &x, const Polynomial_t y) const {
			Polynomial_t a;
			return _pd.areEqual(convert(a, x), y);
		}
		
		bool areEqual(const Polynomial_t y, const Element &x) const {
			Polynomial_t a;
			return _pd.areEqual(convert(a, x), y);
		}
		
		bool areEqual(const Polynomial_t &x, const Polynomial_t y) const {
			return _pd.areEqual(x, y);
		}

		Element &normalizeIn(Element &x) const {
			if (isZero(x)) {
				return x;
			}
			
			_pd.assign(x, one);
			return x;
		}
		
		Element &normalize(Element &z, const Element &x) const {
			_pd.assign(z, x);
			return normalizeIn(z);
		}
		
		bool areAssociates(const Element &x, const Element &y) const {
			bool xz = isZero(x);
			bool yz = isZero(y);
			
			return (xz && yz) || (!xz && !yz && x.e == y.e);
		}

		bool isZero(const Element &x) const {
			return _pd.isZero(x.f);
		}

		bool isOne(const Element &x) const {
			return x.e == 0 && _pd.isOne(x.f);
		}

		bool isMOne(const Element &x) const {
			return x.e == 0 && _pd.isMOne(x.f);
		}

		std::ostream &write(std::ostream &os) const {
			return _pd.write(os << "polynomials mod (", _f) << ")^" << _e;
		}

		std::istream &read(std::istream &is) {
			return is;
		}

		std::ostream &write(std::ostream &os, const Element &x) const {
			Polynomial_t a;
			return _pd.write(os, this->convert(a, x));
		}

		std::ostream &write2(std::ostream &os, const Element &x) const {
			Polynomial_t a;
			_pd.write(os << "(", x.f) << ")*(";
			_pd.write(os, _f) << ")^(" << x.e << ")";
			return os;
		}

		std::istream &read(std::istream &is, Element &x) const {
			Polynomial_t a;
			_pd.read(is, a);
			init(x, a);
			return is;
		}

		Element &addin(Element &x, const Element &y) const {
			Polynomial_t xp, yp;
			
			convert(xp, x);
			convert(yp, y);
			
			_pd.addin(xp, yp);
			
			return init(x, xp);
		}

		Element &add(Element &x, const Element &y, const Element &z) const {
			assign(x, y);
			return addin(x, z);
		}

		Element &subin(Element &x, const Element &y) const {
			Element a;
			neg(a, y);
			return addin(x, a);
		}

		Element &sub(Element &x, const Element &y, const Element &z) const {
			neg(x, z);
			return addin(x, y);
		}

		Element &mulin(Element &x, const Element &y) const {
			if (isZero(y) || x.e + y.e >= _e) {
				return assign(x, zero);
			}
			
			//x.e += y.e;
			//_pd.mulin(x.f, y.f);
			
			Polynomial_t a, b;
			convert(a, x);
			convert(b, y);
			
			_pd.write(std::cout << "a: ", a) << std::endl;
			_pd.write(std::cout << "b: ", b) << std::endl;
			
			_pd.mulin(a, b);
			
			_pd.write(std::cout << "c: ", a) << std::endl;
			
			init(x, a);
			
			write(std::cout << "x: ", x) << std::endl;
			
			return x;
		}

		Element &mul(Element &x, const Element &y, const Element &z) const {
			assign(x, y);
			return mulin(x, z);
		}

		Element &divin(Element &x, const Element &y) const {
			x.e -= y.e;
			
			Polynomial_t g,s,t;
			_pd.gcd(g, s, t, y.f, _fe);
			
			//_pd.write(std::cout << "s: ", s) << std::endl;
			//_pd.write(std::cout << "x: ", x.f) << std::endl;
			
			_pd.mulin(x.f, s);
			
			//_pd.write(std::cout << "x*s: ", x.f) << std::endl;
			
			_pd.modin(x.f, _fe);
			
			//_pd.write(std::cout << "(x*s)%f: ", x.f) << std::endl;
			
			return x;
		}

		Element &div(Element &x, const Element &y, const Element &z) const {
			assign(x, y);
			return divin(x, z);
		}

		Element &negin(Element &x) const {
			_pd.negin(x.f);
			return x;
		}

		Element &neg(Element &x, const Element &y) const {
			x.e = y.e;
			_pd.neg(x.f, y.f);
			return x;
		}

		Element &invin(Element &x) const {
			if (isZero(x) || x.e > 0) {
				return assign(x, zero);
			}
			
			Polynomial_t g,s,t;
			_pd.gcd(g,s,t,x.f,_fe);
			_pd.assign(x.f, s);
			return x;
		}

		Element &inv(Element &x, const Element &y) const {
			assign(x, y);
			return invin(x);
		}

		Element &axpyin(Element &r, const Element &a, const Element &x) const {
			mulin(r, a);
			return addin(r, x);
		}

		Element &axpy(Element &r, const Element &a, const Element &x, const Element &y) const {
			assign(r, a);
			return axpyin(r,x,y);
		}

		// PIR Functions

		bool isDivisor(const Element &a, const Element &b) const {
			return !isZero(b) && a.e >= b.e;
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

		//Element& normalIn(Element &x) const
		//{
		//	return x;
		//}
		
		bool isUnit(const Element& x) const {
			return !_pd.isZero(x.f) && x.e == 0;
		}

		// g = gcd(a,b)
		Element& gcd(Element &g, const Element &a, const Element &b) const {
			if (isZero(a)) {
				return assign(g, b);
			}
			
			if (isZero(b)) {
				return assign(g, a);
			}
						
			_pd.assign(g.f, 1);
			g.e = a.e < b.e ? a.e : b.e;
			
			return g;
		}
		
		Element& gcdin(Element &a, const Element &b) const {
			Element c;
			assign(c, a);
			return gcd(a, c, b);
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

#endif // LINBOX_ring_givaro_poly_local_H
