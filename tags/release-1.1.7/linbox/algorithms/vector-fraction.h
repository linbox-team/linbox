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
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_vector_fraction_H
#define __LINBOX_vector_fraction_H

#include <linbox/linbox-config.h>
#include <stdio.h>
#include <linbox/vector/vector-traits.h>
#undef _D

namespace LinBox 
{
	
	/** utility function to reduce a rational pair to lowest form */
	template<class Domain>
	void reduceIn(Domain& D, std::pair<typename Domain::Element, typename Domain::Element> &frac){
		linbox_check(!D.isZero(frac.second));

		if (D.isZero(frac.first)){
			D.init(frac.second, 1);
			return;
		}

		typename Domain::Element gcd;
		D.gcd(gcd, frac.first, frac.second);
		D.divin(frac.first, gcd);
		D.divin(frac.second, gcd);
	};

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
		const Domain& _D;
		Element zero;

		/** 
		 * constructor from vector of rational numbers
		 * reduces individual pairs in-place first unless alreadyReduced=true
		 */
		VectorFraction(const Domain& D, FVector& frac) // bool alreadyReduced = false)
			: _D(D) {
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
		VectorFraction(const Domain& D, size_t n) 
			: _D(D) {
			D.init(zero, 0);
			D.init(denom, 1);
			numer = Vector(n);
			typename Vector::iterator j;
			
			for (j=numer.begin(); j!=numer.end(); j++)
				D.assign(*j, zero);
		}
		
		/** copy constructor */
		VectorFraction(const VectorFraction<Domain>& VF) 
			: _D(VF._D) {
			copy(VF);
		}

		/** copy without construction */
		void copy(const VectorFraction<Domain>& VF) {
			//assumes _D = VF._D
			denom = VF.denom;
			numer.resize(VF.numer.size());
			typename Vector::iterator i;
			typename Vector::const_iterator j;
		  
			for (i=numer.begin(), j=VF.numer.begin(); i!=numer.end(); i++, j++)
				_D.assign(*i, *j);
		}
		
		/** clear and resize without construction */
		void clearAndResize(size_t size) {
			_D.init(denom, 1);
			typename Vector::iterator i;
			numer.resize(size);
			for (i=numer.begin(); i!=numer.end(); i++)
				_D.init(*i, 0);
		}

		/** 
		 * Replaces *this with a linear combination of *this and other 
		 * such that the result has denominator == gcd(this->denom, other.denom)
		 * see Mulders+Storjohann : 'Certified Dense Linear System Solving' Lemma 2.1
		 * return value of true means that there was some improvement (ie denom was reduced)
		 */
		bool combineSolution(const VectorFraction<Domain>& other) { 
			if (_D.isDivisor(other.denom, denom)) return false;
			if (_D.isDivisor(denom, other.denom)) {
				denom = other.denom;
				numer = other.numer;
				return true;
			}
			Element s, t, g;
			_D.xgcd(g, s, t, denom, other.denom);
			if (_D.areEqual(g, denom)) ; //do nothing
			else {
				denom = g;
				typename Vector::iterator it=numer.begin();
				typename Vector::const_iterator io=other.numer.begin();
				for (; it != numer.end(); it++, io++) {
					_D.mulin(*it, s);
					_D.axpyin(*it, t, *io);
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
		bool boundedCombineSolution(const VectorFraction<Domain>& other, const Element& denBound, Element& g) {

			//this means that new solution won't reduce g
			if (_D.isDivisor(other.denom, g)) return false;
			
			//short-circuit in case the new solution is completely better than old one
			Element _dtmp;
			if (_D.isDivisor(g, _D.gcd(_dtmp, denBound, other.denom))) {
				denom = other.denom;
				numer = other.numer;
				g = _dtmp;
				return true;
			}

			Element A, g2, lincomb;
			_D.gcd(g, other.denom, g); //we know this reduces g

			// find A s.t. gcd(denBound, denom + A*other.denom) = g
			// strategy: pick random values of A <= d(y_0)
			integer tmp;
			_D.convert(tmp, denBound);
			typename Domain::RandIter randiter(_D, tmp); //seed omitted
			// TODO: I don't think this random iterator has high-quality low order bits, which are needed
			do {
				randiter.random(A);
				_D.assign(lincomb, denom);
				_D.axpyin(lincomb, A, other.denom);
				_D.gcd(g2, lincomb, denBound);
			}
			while (!_D.areEqual(g, g2));
			
			_D.assign(denom, lincomb);
			typename Vector::iterator it=numer.begin();
			typename Vector::const_iterator io=other.numer.begin();
			for (; it != numer.end(); it++, io++) 
				_D.axpyin(*it, A, *io);
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
					const Element& n2, const Element d2) {
			//this means that new solution won't reduce g
			if (_D.isDivisor(d1, d2)) return false;
			
			//short-circuit in case the new solution is completely better than old one
			if (_D.isDivisor(d2, d1)) {
				copy(other);
				n1 = n2;
				d1 = d2;
				return true;
			}

			Element A, g, l, n1d2_g, n2d1_g, lincomb, g2, tmpe, one;

			_D.gcd(g, d1, d2);   //compute gcd
			_D.mul(l, d1, d2);
			_D.divin(l, g);      //compute lcm
			
			_D.div(n1d2_g, d2, g);
			_D.mulin(n1d2_g, n1);   //compute n1.d2/g
			_D.div(n2d1_g, d1, g);
			_D.mulin(n2d1_g, n2);   //compute n2.d1/g

			// find A s.t. gcd(denBound, denom + A*other.denom) = g
			// strategy: pick random values of A <= lcm(d(denom), d(other.denom))
			integer tmp;
			_D.mul(tmpe, denom, other.denom);
			_D.convert(tmp, tmpe);
			_D.init(one, 1);
			typename Domain::RandIter randiter(_D, tmp); //seed omitted
			// TODO: I don't think this random iterator has high-quality low order bits, which are needed
			do {
				randiter.random(A);
				_D.assign(lincomb, n1d2_g);
				_D.axpyin(lincomb, A, n2d1_g);
				_D.gcd(g2, lincomb, l);
			}
			while (!_D.areEqual(one, g2));
			
			this->axpyin(A, other);
			_D.lcmin(d1, d2);

			return true;
		}

		/**
		 * this += a * x.   performs a rational axpy with an integer multiplier
		 * returns (*this)
		 */
		VectorFraction<Domain>& axpyin(Element& a, const VectorFraction<Domain>& x) {
			Element a_prime, gcd_a_xdenom, xdenom_prime;
			_D.gcd(gcd_a_xdenom, a, x.denom);
			_D.div(a_prime, a, gcd_a_xdenom);
			_D.div(xdenom_prime, x.denom, gcd_a_xdenom);
			
			Element cdf; //common denominator factor; multiply both sides by this and divide at end
			_D.gcd(cdf, denom, xdenom_prime);
			_D.divin(denom, cdf);
			_D.divin(xdenom_prime, cdf);

			// we perform numer[i] = xdenom_prime * numer[i] + a_prime * denom * x.denom[i]
			// so multiply denom into a_prime and save a multiplication on each entry
			_D.mulin(a_prime, denom);

			typename Vector::iterator i = this->numer.begin();
			typename Vector::const_iterator j = x.numer.begin();
			for (; i != this->numer.end(); i++, j++) {
				_D.mulin(*i, xdenom_prime);
				_D.axpyin(*i, a_prime, *j);
			}
			
			_D.mulin(denom, cdf);
			_D.mulin(denom, xdenom_prime);
			simplify();
			return *this;
		}

		/** write to a stream */
		std::ostream& write(std::ostream& os) const {
			os << "[";
			for (typename Vector::const_iterator it=numer.begin(); it != numer.end(); it++) {
				if (it != numer.begin()) os << " ";
				os << *it;
			}
			return os << "]/" << denom;
		}
		
		/** convert to 'answer' type of lifting container */
		FVector& toFVector(FVector& result) const {
			linbox_check(numer.size()==result.size());
			typename Vector::const_iterator it=numer.begin();
			typename FVector::iterator ir=result.begin();
			for (; it != numer.end(); it++, ir++) {
				_D.assign(ir->first, *it);
				_D.assign(ir->second, denom);
			}
			return result;
		}

		/** reduces to simplest form, returns (*this) */
		VectorFraction<Domain>& simplify() { 
			typename Vector::iterator i;
			Element gcd;
			_D.init(gcd, denom);
			vectorGcdIn(gcd, _D, numer);
			
			_D.divin(denom, gcd);
			for (i=numer.begin(); i!=numer.end(); i++) 
				_D.divin(*i, gcd);
			return (*this);
		}
	};

}

#endif //__LINBOX_vector_fraction_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
