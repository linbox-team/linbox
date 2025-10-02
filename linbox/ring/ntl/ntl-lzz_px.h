/* Copyright (C) 2005 LinBox
 * Copyright (C) 2011 LinBox
 *
 *
 * Written by Daniel Roche, August 2005
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

/*! @file ring/ntl/ntl-lzz_pX.h
 * @ingroup ring
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_lzz_px_H
#define __LINBOX_field_ntl_lzz_px_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <vector>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pXFactoring.h>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include <givaro/zring.h>
#include "linbox/ring/ntl/ntl-lzz_p.h"
#include "linbox/ring/ntl/ntl-lzz_pe.h"
#include "linbox/integer.h"


// Namespace in which all LinBox code resides
namespace LinBox
{
	class NTL_zz_pX_Initialiser {
	public :
		NTL_zz_pX_Initialiser( const Integer & q, size_t e = 1) {
			linbox_check(e == 1);
			if ( q > 0 )
				NTL::zz_p::init(int64_t(q)); // it's an error if q not prime, e not 1
		}

		// template <class ElementInt>
		// NTL_zz_pX_Initialiser(const ElementInt& d) {
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		NTL_zz_pX_Initialiser () { }

	};

	/** Ring (in fact, a unique factorization domain) of polynomial with
	 * coefficients in class NTL_zz_p (integers mod a wordsize prime).
	 * All the same functions as any other ring, with the addition of:
	 * Coeff (type), CoeffField (type), getCoeffField, setCoeff, getCoeff,
	 * leadCoeff, deg
	 */
	class NTL_zz_pX :  public NTL_zz_pX_Initialiser, public Givaro::UnparametricOperations<NTL::zz_pX> {
	public:
		typedef NTL::zz_pX Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		typedef NTL_zz_pE QuotientRing;

		typedef NTL_zz_p CoeffField;
		typedef NTL::zz_p Coeff;
		// typedef NTL::zz_pX Element;

		const Element zero,one,mOne ;


		/** Standard LinBox field constructor.  The paramters here
		 * (prime, exponent) are only used to initialize the coefficient field.
		 */
		NTL_zz_pX( const integer& p, size_t e = 1 ) :
			// Givaro::ZRing<NTL::zz_pX>(p, e), _CField(p,e)
			NTL_zz_pX_Initialiser(p,e),Father_t ()
			, zero( NTL::to_zz_pX(0)),one( NTL::to_zz_pX(1)),mOne(-one)
			, _CField(p,e)
		{}

		/** Constructor from a coefficient field */
		NTL_zz_pX( CoeffField cf ) :
			NTL_zz_pX_Initialiser(cf.cardinality()),Father_t ()
			,zero( NTL::to_zz_pX(0)),one( NTL::to_zz_pX(1)),mOne(-one)
			,_CField(cf)
		{}

		Element& init( Element& p) const
		{	return init(p, 0); }

		Element& init( Element& p, integer n ) const
		{
			p = 0;
			integer base;
			_CField.cardinality(base);
			for (int i = 0; n > 0; n /= base, ++i)
				NTL::SetCoeff( p, i, int64_t(n%base) );
			return p;
		}

		integer& convert( integer& n, Element& p ) const
		{
			integer base;
			Coeff x; integer a;
			_CField.cardinality(base);
			n = 0;
			for (int i = (int)deg(p); i >= 0; --i)
			{
				n *= base;
				NTL::GetCoeff(x, p, i );
				n += _CField.convert(a, x);
			}
			return n;
		}

#if 0
		/** Initialize p to the constant y (p = y*x^0) */
		template <class ANY>
		Element& init( Element& p, const ANY& y = 0) const
		{
			Coeff temp;
			_CField.init( temp, y );
			return p = temp;
		}
#endif

		/** Initialize p to the constant y (p = y*x^0) */
		Element& init( Element& p, const Coeff& y ) const
		{
			return p = y;
		}

		/** Initialize p from a vector of coefficients.
		 * The vector should be ordered the same way NTL does it: the front
		 * of the vector corresponds to the trailing coefficients, and the back
		 * of the vector corresponds to the leading coefficients.  That is,
		 * v[i] = coefficient of x^i.
		 */
		template <class ANY1, class ANY2, template <class T1, class T2> class Vect>
		Element& init( Element& p, const Vect<ANY1,ANY2>& v ) const
		{
			p = 0;
			Coeff temp;
			for( long i = 0; i < (long)v.size(); ++i ) {
				_CField.init( temp, v[ (size_t) i ] );
				if( !_CField.isZero(temp) )
					NTL::SetCoeff( p, i, temp );
			}
			return p;
		}

		template <class ANY>
		Element& init( Element& p, const BlasSubvector<ANY>& v ) const
		{
			p = 0;
			Coeff temp;
			for( long i = 0; i < (long)v.size(); ++i ) {
				_CField.init( temp, v[ (size_t) i ] );
				if( !_CField.isZero(temp) )
					NTL::SetCoeff( p, i, temp );
			}
			return p;
		}


		/** Initialize p from a vector of coefficients.
		 * The vector should be ordered the same way NTL does it: the front
		 * of the vector corresponds to the trailing coefficients, and the back
		 * of the vector corresponds to the leading coefficients.  That is,
		 * v[i] = coefficient of x^i.
		 */
		Element& init( Element& p, const std::vector<Coeff>& v ) const
		{
			p = 0;
			for( long i = 0; i < (long)v.size(); ++i )
                NTL::SetCoeff( p, i, v[ (size_t) i ] );
			return p;
		}

		/** Convert p to a vector of coefficients.
		 * The vector will be ordered the same way NTL does it: the front
		 * of the vector corresponds to the trailing coefficients, and the back
		 * of the vector corresponds to the leading coefficients.  That is,
		 * v[i] = coefficient of x^i.
		 */
		template< class ANY >
		std::vector<ANY>& convert( std::vector<ANY>& v, const Element& p ) const
		{
			v.clear();
			ANY temp;
			for( long i = 0; i <= (long)this->deg(p); ++i ) {
				_CField.convert( temp, NTL::coeff( p, i ) );
				v.push_back( temp );
			}
			return v;
		}

		/** Convert p to a vector of coefficients.
		 * The vector will be ordered the same way NTL does it: the front
		 * of the vector corresponds to the trailing coefficients, and the back
		 * of the vector corresponds to the leading coefficients.  That is,
		 * v[i] = coefficient of x^i.
		 */
		std::vector<Coeff>& convert( std::vector<Coeff>& v, const Element& p ) const
		{
			v.clear();
			for( long i = 0; i <= (long)this->deg(p); ++i )
				v.push_back( NTL::coeff(p,i) );
			return v;
		}

		/** Test if an element equals zero */
		bool isZero( const Element& x ) const
		{
			return NTL::IsZero(x);
		}

		/** Test if an element equals one */
		bool isOne( const Element& x ) const
		{
			return NTL::IsOne(x);
		}

		bool isUnit( const Element& x ) const
		{
			return ( (this->deg(x) == 0) &&
				 ( _CField.isUnit( NTL::ConstTerm(x) ) ) );
		}

        bool isMOne(const Element& x) const
		{
			return ( (this->deg(x) == 0) &&
				 ( _CField.isMOne( NTL::ConstTerm(x) ) ) );

		}

		bool isIrreducible(const Element &x) const {
			return NTL::DetIrredTest(x);
		}

		void factor(std::vector<std::pair<Element, long>> &factors, const Element &f) const {
			NTL::Vec<NTL::Pair<Element, long>> factors1;
			NTL::CanZass(factors1, f);

			for (int i = 1; i <= factors1.length(); i++) {
				NTL::Pair<Element, long> tup = factors1(i);
				std::pair<Element, long> tmp(tup.a, tup.b);
				factors.push_back(tmp);
			}
		}

		void squareFree(std::vector<std::pair<Element, long>> &factors, const Element &f) const {
			NTL::Vec<NTL::Pair<Element, long>> factors1;
			NTL::SquareFreeDecomp(factors1, f);

			for (int i = 1; i <= factors1.length(); i++) {
				NTL::Pair<Element, long> tup = factors1(i);
				std::pair<Element, long> tmp(tup.a, tup.b);
				factors.push_back(tmp);
			}
		}

		/** The LinBox field for coefficients */
		const CoeffField& getCoeffField() const
		{ return _CField; }

		/** Get the degree of a polynomial
		 * Unlike NTL, deg(0)=0.
		 */
		size_t deg( const Element& p ) const
		{
			long temp = NTL::deg(p);
			if( temp == -1 ) return 0;
			else return static_cast<size_t>(temp);
		}

		/** r will be set to the reverse of p. */
		Element& rev( Element& r, const Element& p ) {
			NTL::reverse(r,p);
			return r;
		}

		/** r is itself reversed. */
		Element& revin( Element& r ) {
			return r = NTL::reverse(r);
		}

		/** Get the leading coefficient of this polynomial. */
		Coeff& leadCoeff( Coeff& c, const Element& p ) const
		{
			c = NTL::LeadCoeff(p);
			return c;
		}

		Element& monic(Element& r, const Element& p) const {
			r = p;
			NTL::MakeMonic(r);
			return r;
		}

		Element& monicIn(Element& p) const {
			NTL::MakeMonic(p);
			return p;
		}

		/** Get the coefficient of x^i in a given polynomial */
		Coeff& getCoeff( Coeff& c, const Element& p, size_t i ) const
		{
			c = NTL::coeff( p, (long)i );
			return c;
		}

		/** Set the coefficient of x^i in a given polynomial */
		Element& setCoeff( Element& p, size_t i, const Coeff& c ) const
		{
			NTL::SetCoeff(p,(long)i,c);
			return p;
		}

		Element &leftShift(Element &x, const Element &a, size_t shift) const {
			x = a << shift;
			return x;
		}

		Element &leftShiftIn(Element &a, size_t shift) const {
			a <<= shift;
			return a;
		}

		Element &rightShift(Element &x, const Element &a, size_t shift) const {
			x = a >> shift;
			return x;
		}

		Element &rightShiftIn(Element &a, size_t shift) const {
			a >>= shift;
			return a;
		}

		Element& mulCoeffIn(Element &p, const Coeff &c) const {
			p *= c;
			return p;
		}

		Element& pow(Element& x, const Element& a, long e) const {
			NTL::power(x, a, e);
			return x;
		}

		/** Get the quotient of two polynomials */
		Element& quo( Element& res, const Element& a, const Element& b ) const
		{
			NTL::div(res,a,b);
			return res;
		}

		/** a = quotient of a, b */
		Element& quoin( Element& a, const Element& b ) const
		{
			return a /= b;
		}

		/** Get the remainder under polynomial division */
		Element& rem( Element& res, const Element& a, const Element& b ) const
		{
			NTL::rem(res,a,b);
			return res;
		}

		/** a = remainder of a,b */
		Element& remin( Element& a, const Element& b ) const
		{
			return a %= b;
		}

		/** Get the quotient and remainder under polynomial division */
		void quorem( Element& q, Element& r,
			     const Element& a, const Element& b ) const
		{
			NTL::DivRem(q,r,a,b);
		}

		bool isDivisor(const Element &a, const Element &b) const {
			if (isZero(b)) {
				return false;
			}

			Element tmp;
			rem(tmp, a, b);
			return isZero(tmp);
		}

		// a = b^(-1) % f
		Element& invMod(Element &a, const Element &b, const Element &f) const
		{
			NTL::InvMod(a, b, f);
			return a;
		}

		Element& inv( Element& y, const Element& x ) const
		{
			// Element one(0, 1);
			return quo(y,one,x);
		}

		Element& invin( Element& y ) const
		{
			Element x = y;
			return inv(y, x);
		}

		/** Get the greatest commonn divisor of two polynomials */
		Element& gcdin( Element& a, const Element& b ) const
		{
			Element res;
			NTL::GCD(res,a,b);
			a = res;
			return a;
		}

		Element& gcd( Element& res, const Element& a, const Element& b ) const
		{
			NTL::GCD(res,a,b);
			return res;
		}

		Element& gcd( Element& res, Element& s, Element& t, const Element& a, const Element& b ) const
		{
			NTL::XGCD(res,s,t,a,b);
			return res;
		}

		Element &dxgcd(Element &g, Element &s, Element &t, Element &u, Element &v, const Element &a, const Element &b) const {
			gcd(g,s,t,a,b);

			div(u,a,g);
			div(v,b,g);

			return g;
		}

		/** Get the least common multiple of two polynomials */
		Element& lcmin( Element& a, const Element& b ) const
		{
			Element tmp, res;
			gcd(tmp, a, b);
			div(res, a, tmp);
			return mul(a, res, b);
		}

		Element& lcm( Element& res, const Element& a, const Element& b ) const
		{
            if (isZero(a) || isZero(b)) return assign(res,zero);
			Element tmp;
			gcd(tmp,a,b);
			div(res, a, tmp);
			return mulin(res, b);
		}

		/** Get characteristic of the field - same as characteristic of
		 * coefficient field. */
		integer& characteristic( integer& c ) const
		{ return _CField.characteristic(c); }

		/** Get the cardinality of the field.  Since the cardinality is
		 * infinite, by convention we return -1.
		 */
		integer& cardinality( integer& c ) const
		{ return c = static_cast<integer>(-1); }

		static inline integer maxCardinality()
		{ return CoeffField::maxCardinality(); }
		/** Write a description of the field */
		// Oustide of class definition so write(ostream&,const Element&) from
		// Givaro::ZRing still works.

		std::ostream& write( std::ostream& os ) const
		{
			return os << "Polynomial ring using NTL::zz_pX";
		}

		std::ostream& write( std::ostream& os, const Element& x) const {
			// return Father_t::write(os, x);
			if (isZero(x)) {
				os << "0";
				return os;
			}

			bool first = true;
			for (size_t i = 0; i <= deg(x); i++) {
				Coeff xi = NTL::coeff(x, i);
				if (xi != 0) {
					if (!first) {
						os << "+";
					}

					if (xi == 1 && i > 0) {
						os << "x";
					} else {
						os << xi;
						if (i > 0) {
							os << "*x";
						}
					}

					if (i > 1) {
						os << "^" << i;
					}
					first = false;
				}
			}
			return os;
		}

		std::istream& read(std::istream& i, Element& p) const {
			long deg;
			i >> deg;

			std::vector<integer> coeffs(deg + 1);
			for(; deg >= 0; --deg) {
				i >> coeffs[deg];
			}

			init(p, coeffs);
			return i;
		}

		/** Conversion to scalar types doesn't make sense and should not be
		 * used.  Use getCoeff or leadCoeff to get the scalar values of
		 * specific coefficients, and then convert them using coeffField()
		 * if needed.
		 */
		template< class ANY >
		ANY& convert( ANY& x, const Element& y ) const
		{ return x; }

        const NTL_zz_pE quotient(const Element &f) const {
	integer c;
	QuotientRing QR(characteristic(c), f);
	return QR;
        }

	protected:
		CoeffField _CField;
	}; // end of class NTL_zz_pX




	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::zz_pX> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::zz_pX> {
	public:
		typedef NTL::zz_pX Element;
        typedef size_t Residu_t;

		UnparametricRandIter(const NTL_zz_pX & F ,
                             const uint64_t& seed = 0,
                             const Residu_t& size = 0) :
            _size(size), _seed(seed), _ring(F) {
			if(_seed == 0)
                NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(std::time(nullptr))));
			else
                NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(_seed)));
		}

        const NTL_zz_pX& ring() const { return _ring; }
		UnparametricRandIter(const UnparametricRandIter<NTL::zz_pX>& R) :
            _size(R._size), _seed(R._seed), _ring(R._ring) {
			if(_seed == 0)
                NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(std::time(nullptr))));
            else
                NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(_seed)));
		}

		Element& random (Element& x) const {
			NTL::random(x, 1);
			return x;
		}

		Element& random(Element &x, size_t d) const {
			NTL::random(x, d + 1);
			return x;
		}

		Element& randomIrreducible(Element &x, size_t d) const {
			NTL::BuildIrred(x, (long) d);
			return x;
		}

	protected:
		size_t _size;
		uint64_t _seed;
        const NTL_zz_pX& _ring;
	}; // class UnparametricRandIters

	template<>
	class Hom<NTL_zz_pX, NTL_zz_pE> {
	public:
		typedef NTL_zz_pX Source;
		typedef NTL_zz_pE Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;

		Hom(const Source& S, const Target& T) : _source(S), _target(T) {}

		Elt& image(Elt& t, const SrcElt& s) {
			return t = NTL::conv<NTL::zz_pE>(s);
		}

		SrcElt& preimage(SrcElt& s, const Elt& t) {
			return s = NTL::conv<NTL::zz_pX>(t);
		}

		const Source& source() { return _source;}
		const Target& target() { return _target;}
	private:
		const Source& _source;
		const Target& _target;
	}; // end Hom<NTL_zz_pX, NTL_zz_pE>

} // end of namespace LinBox

#endif // __LINBOX_field_ntl_lzz_px_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
