/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/vector-fraction.h
 * Copyright (C) 2004 David Pritchard
 *
 * Written by David Pritchard <daveagp@mit.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_vector_fraction_H
#define __LINBOX_vector_fraction_H

#include "linbox/linbox-config.h"
#include <stdio.h>
#include "linbox/vector/vector-traits.h"

namespace LinBox
{

	/** utility function to reduce a rational pair to lowest form */
	template<class Domain>
	void reduceIn(Domain& D, std::pair<typename Domain::Element, typename Domain::Element> &frac)
	{
		linbox_check(!D.isZero(frac.second));

		if (D.isZero(frac.first)){
			D.init(frac.second, 1);
			return;
		}

		typename Domain::Element gcd;
		D.gcd(gcd, frac.first, frac.second);
		D.divin(frac.first, gcd);
		D.divin(frac.second, gcd);
	}

	/** utility function to gcd-in a vector of elements over a domain */
	//this could be replaced by a fancier version that combines elements linearly at random
	template<class Domain, class Vector>
	void vectorGcdIn(typename Domain::Element& result, Domain& D, Vector& v) {
		for (typename Vector::iterator i = v.begin(); i != v.end(); i++)
			D.gcdin(result, *i);
	}

	/** utility function, returns gcd of a vector of elements over a domain */
	// this could be replaced by a fancier version that combines elements linearly at random
	template<class Domain, class Vector>
	typename Domain::Element vectorGcd(Domain& D, Vector& v) {
		typename Domain::Element result;
		D.init(result, 0);
		vectorGcdIn(result, D, v);
		return result;
	}


	/**
	 * \brief VectorFraction<Domain> is a vector of rational elements with common reduced denominator.
	 * Here Domain is a ring supporting the gcd, eg NTL_ZZ or PID_integer
	 *	 For compatability with the return type of rationalSolver, it allows conversion from/to
	 *       std::vector<std::pair<Domain::Element> >.
	 *       All functions will return the fraction in reduced form, calling reduce() if necessary.
	 */

	template<class Domain>
	class VectorFraction{
	public:
		typedef typename Domain::Element Element;
		typedef typename std::pair<Element, Element> Fraction;
		typedef typename std::vector<Fraction> FVector;
		typedef typename Vector<Domain>::Dense Vector;

		Vector numer;
		Element denom;
		const Domain& _domain;
		Element zero;

		/**
		 * constructor from vector of rational numbers
		 * reduces individual pairs in-place first unless alreadyReduced=true
		 */
		VectorFraction(const Domain& D, FVector& frac
			       //,bool alreadyReduced = false
			      ) :
			_domain(D)
		{
			bool alreadyReduced = false;
			typename FVector::iterator i;

			D.init(zero, 0);
			D.init(denom, 1);
			if (!alreadyReduced)
				for (i=frac.begin(); i!=frac.end(); i++)
					reduceIn(D, *i);

			for (i=frac.begin(); i!=frac.end(); i++) {
				linbox_check(!D.isZero(i->second));
				D.lcmin(denom, i->second);
			}

			numer = Vector(frac.size());
			typename Vector::iterator j;

			for (i=frac.begin(), j=numer.begin(); i!=frac.end(); i++, j++){
				D.mul(*j, denom, i->first);
				D.divin(*j, i->second);
			}
		}

		/** allocating constructor, returns [0, 0, ... 0]/1 */
		VectorFraction(const Domain& D, size_t n) :
			_domain(D)
		{
			D.init(zero, 0);
			D.init(denom, 1);
			numer = Vector(n);
			typename Vector::iterator j;

			for (j=numer.begin(); j!=numer.end(); j++)
				D.assign(*j, zero);
		}

		/** copy constructor */
		VectorFraction(const VectorFraction<Domain>& VF) :
			_domain(VF._domain)
		{
			copy(VF);
		}

		/** copy without construction */
		void copy(const VectorFraction<Domain>& VF)
		{
			//assumes _domain = VF._domain
			denom = VF.denom;
			numer.resize(VF.numer.size());
			typename Vector::iterator i;
			typename Vector::const_iterator j;

			for (i=numer.begin(), j=VF.numer.begin(); i!=numer.end(); i++, j++)
				_domain.assign(*i, *j);
		}

		/** clear and resize without construction */
		void clearAndResize(size_t size)
		{
			_domain.init(denom, 1);
			typename Vector::iterator i;
			numer.resize(size);
			for (i=numer.begin(); i!=numer.end(); i++)
				_domain.init(*i, 0);
		}

		/**
		 * Replaces *this with a linear combination of *this and other
		 * such that the result has denominator == gcd(this->denom, other.denom)
		 * see Mulders+Storjohann : 'Certified Dense Linear System Solving' Lemma 2.1
		 * return value of true means that there was some improvement (ie denom was reduced)
		 */
		bool combineSolution(const VectorFraction<Domain>& other)
		{
			if (_domain.isDivisor(other.denom, denom)) return false;
			if (_domain.isDivisor(denom, other.denom)) {
				denom = other.denom;
				numer = other.numer;
				return true;
			}
			Element s, t, g;
			_domain.xgcd(g, s, t, denom, other.denom);
			if (_domain.areEqual(g, denom)) ; //do nothing
			else {
				denom = g;
				typename Vector::iterator it=numer.begin();
				typename Vector::const_iterator io=other.numer.begin();
				for (; it != numer.end(); it++, io++) {
					_domain.mulin(*it, s);
					_domain.axpyin(*it, t, *io);
				}
				return true;
			}
			return false;
		}

		/**
		 * Adds in-place to *this a multiple of other
		 * such that the result has gcd(denominator, denBound) == gcd(this->denom, other.denom, denBound)
		 * see Mulders+Storjohann : 'Certified Dense Linear System Solving' Lemma 6.1
		 * return value of true means that there was some improvement (ie gcd(denom, denBound) was reduced)
		 * g is gcd(denom, denBound), and is updated by this function when there is improvement
		 */
		bool boundedCombineSolution(const VectorFraction<Domain>& other, const Element& denBound, Element& g)
		{

			//this means that new solution won't reduce g
			if (_domain.isDivisor(other.denom, g)) return false;

			//short-circuit in case the new solution is completely better than old one
			Element _dtmp;
			if (_domain.isDivisor(g, _domain.gcd(_dtmp, denBound, other.denom))) {
				denom = other.denom;
				numer = other.numer;
				g = _dtmp;
				return true;
			}

			Element A, g2, lincomb;
			_domain.gcd(g, other.denom, g); //we know this reduces g

			// find A s.t. gcd(denBound, denom + A*other.denom) = g
			// strategy: pick random values of A <= d(y_0)
			integer tmp;
			_domain.convert(tmp, denBound);
			typename Domain::RandIter randiter(_domain, tmp); //seed omitted
			// TODO: I don't think this random iterator has high-quality low order bits, which are needed
			do {
				randiter.random(A);
				_domain.assign(lincomb, denom);
				_domain.axpyin(lincomb, A, other.denom);
				_domain.gcd(g2, lincomb, denBound);
			}
			while (!_domain.areEqual(g, g2));

			_domain.assign(denom, lincomb);
			typename Vector::iterator it=numer.begin();
			typename Vector::const_iterator io=other.numer.begin();
			for (; it != numer.end(); it++, io++)
				_domain.axpyin(*it, A, *io);
			return true;
		}

		/**
		 * Adds in-place to *this a multiple of other to create an improved certificate ("z")
		 * n1/d1 = *this . b, n2/d2 = other . b   in reduced form
		 * n1/d1 are updated so that new denominator is lcm(d1, d2);
		 * see Mulders+Storjohann : 'Certified Dense Linear System Solving' Lemma 6.2
		 * return value of true means that there was some improvement (ie d1 was increased)
		 */
		bool combineCertificate(const VectorFraction<Domain>& other, Element& n1, Element& d1,
					const Element& n2, const Element d2)
		{
			//this means that new solution won't reduce g
			if (_domain.isDivisor(d1, d2)) return false;

			//short-circuit in case the new solution is completely better than old one
			if (_domain.isDivisor(d2, d1)) {
				copy(other);
				n1 = n2;
				d1 = d2;
				return true;
			}

			Element A, g, l, n1d2_g, n2d1_g, lincomb, g2, tmpe, one;

			_domain.gcd(g, d1, d2);   //compute gcd
			_domain.mul(l, d1, d2);
			_domain.divin(l, g);      //compute lcm

			_domain.div(n1d2_g, d2, g);
			_domain.mulin(n1d2_g, n1);   //compute n1.d2/g
			_domain.div(n2d1_g, d1, g);
			_domain.mulin(n2d1_g, n2);   //compute n2.d1/g

			// find A s.t. gcd(denBound, denom + A*other.denom) = g
			// strategy: pick random values of A <= lcm(d(denom), d(other.denom))
			integer tmp;
			_domain.mul(tmpe, denom, other.denom);
			_domain.convert(tmp, tmpe);
			_domain.init(one, 1);
			typename Domain::RandIter randiter(_domain, tmp); //seed omitted
			// TODO: I don't think this random iterator has high-quality low order bits, which are needed
			do {
				randiter.random(A);
				_domain.assign(lincomb, n1d2_g);
				_domain.axpyin(lincomb, A, n2d1_g);
				_domain.gcd(g2, lincomb, l);
			}
			while (!_domain.areEqual(one, g2));

			this->axpyin(A, other);
			_domain.lcmin(d1, d2);

			return true;
		}

		/**
		 * this += a * x.   performs a rational axpy with an integer multiplier
		 * returns (*this)
		 */
		VectorFraction<Domain>& axpyin(Element& a, const VectorFraction<Domain>& x)
		{
			Element a_prime, gcd_a_xdenom, xdenom_prime;
			_domain.gcd(gcd_a_xdenom, a, x.denom);
			_domain.div(a_prime, a, gcd_a_xdenom);
			_domain.div(xdenom_prime, x.denom, gcd_a_xdenom);

			Element cdf; //common denominator factor; multiply both sides by this and divide at end
			_domain.gcd(cdf, denom, xdenom_prime);
			_domain.divin(denom, cdf);
			_domain.divin(xdenom_prime, cdf);

			// we perform numer[i] = xdenom_prime * numer[i] + a_prime * denom * x.denom[i]
			// so multiply denom into a_prime and save a multiplication on each entry
			_domain.mulin(a_prime, denom);

			typename Vector::iterator i = this->numer.begin();
			typename Vector::const_iterator j = x.numer.begin();
			for (; i != this->numer.end(); i++, j++) {
				_domain.mulin(*i, xdenom_prime);
				_domain.axpyin(*i, a_prime, *j);
			}

			_domain.mulin(denom, cdf);
			_domain.mulin(denom, xdenom_prime);
			simplify();
			return *this;
		}

		/** write to a stream */
		std::ostream& write(std::ostream& os) const
		{
			os << "[";
			for (typename Vector::const_iterator it=numer.begin(); it != numer.end(); it++) {
				if (it != numer.begin()) os << " ";
				os << *it;
			}
			return os << "]/" << denom;
		}

		/** convert to 'answer' type of lifting container */
		FVector& toFVector(FVector& result) const
		{
			linbox_check(numer.size()==result.size());
			typename Vector::const_iterator it=numer.begin();
			typename FVector::iterator ir=result.begin();
			for (; it != numer.end(); it++, ir++) {
				_domain.assign(ir->first, *it);
				_domain.assign(ir->second, denom);
			}
			return result;
		}

		/** reduces to simplest form, returns (*this) */
		VectorFraction<Domain>& simplify()
		{
			typename Vector::iterator i;
			Element gcd;
			_domain.init(gcd, denom);
			vectorGcdIn(gcd, _domain, numer);

			_domain.divin(denom, gcd);
			for (i=numer.begin(); i!=numer.end(); i++)
				_domain.divin(*i, gcd);
			return (*this);
		}
	};

}

#endif //__LINBOX_vector_fraction_H
