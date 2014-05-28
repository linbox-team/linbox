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

/*! @file field/NTL/ntl-lzz_pEX.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_lzz_pEX_H
#define __LINBOX_field_ntl_lzz_pEX_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <vector>
#include <NTL/lzz_pEX.h>


#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/field/unparametric.h"
#include "linbox/field/NTL/ntl-lzz_pE.h"
#include "linbox/integer.h"


// Namespace in which all LinBox code resides
namespace LinBox
{
	class NTL_zz_pEX_Initialiser {
	public :
		NTL_zz_pEX_Initialiser( const Integer & q, size_t e = 1) {
			if ( q > 0 )
				NTL::zz_p::init(q); // it's an error if q not prime
			NTL::zz_pX irredPoly = NTL::BuildIrred_zz_pX ((long) e);
			NTL::zz_pE::init(irredPoly); 
		}

		// template <class ElementInt>
		// NTL_zz_pEX_Initialiser(const ElementInt& d) {
			// NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		NTL_zz_pEX_Initialiser () { }

	};

	/** Ring (in fact, a unique factorization domain) of polynomial with
	 * coefficients in class NTL_zz_p (integers mod a wordsize prime).
	 * All the same functions as any other ring, with the addition of:
	 * Coeff (type), CoeffField (type), getCoeffField, setCoeff, getCoeff,
	 * leadCoeff, deg
	 */
	class NTL_zz_pEX :  public NTL_zz_pEX_Initialiser, public FFPACK::UnparametricOperations<NTL::zz_pEX> {
	public:
		typedef NTL::zz_pEX Element ;
		typedef FFPACK::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;


		typedef NTL_zz_pE CoeffField;
		typedef NTL::zz_pE Coeff;
		// typedef NTL::zz_pEX Element;

		const Element zero,one,mOne ;


		/** Standard LinBox field constructor.  The paramters here
		 * (prime, exponent) are only used to initialize the coefficient field.
		 */
		NTL_zz_pEX( const integer& p, size_t e = 1 ) :
			// UnparametricField<NTL::zz_pEX>(p, e), _CField(p,e)
			NTL_zz_pEX_Initialiser(p,e),Father_t ()
			, zero( NTL::to_zz_pEX(0)),one( NTL::to_zz_pEX(1)),mOne(-one)
			, _CField(p,e)
		{}

		/** Constructor from a coefficient field */
		NTL_zz_pEX( CoeffField cf ) :
			NTL_zz_pEX_Initialiser(cf.cardinality()),Father_t ()
			,zero( NTL::to_zz_pEX(0)),one( NTL::to_zz_pEX(1)),mOne(-one)
			,_CField(cf)
		{}

		Element& init( Element& p) const
		{	return init(p, 0); }

		Element& init( Element& p, integer n ) const
		{
			p = 0;
			integer base;
			Coeff a;
			_CField.cardinality(base);
			for (int i = 0; n > 0; n /= base, ++i)
			{
			    _CField.init(a, n%base);
				NTL::SetCoeff( p, i, a );
			}
			return p;
		}

		integer& convert( integer& n, Element& p ) const
		{
			integer base;
			integer d;
			_CField.cardinality(base);
			n = 0;
			for (int i = deg(p); i >= 0; --i)
			{
				n *= base;
				n += _CField.convert(d, NTL::coeff(p, i)); 
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
		template <class ANY>
		Element& init( Element& p, const std::vector<ANY>& v ) const
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
		Element& init( Element& p, const BlasVector<ANY>& v ) const
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
			return ( (this->deg(x) == 0) &&
				 ( _CField.isZero( NTL::ConstTerm(x) ) ) );
		}

		/** Test if an element equals one */
		bool isOne( const Element& x ) const
		{
			return ( (this->deg(x) == 0) &&
				 ( _CField.isOne( NTL::ConstTerm(x) ) ) );
		}

	bool isMOne (const Element& x) const
		{
			return ( (this->deg(x) == 0) &&
				 ( _CField.isMOne( NTL::ConstTerm(x) ) ) );

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

		Element& inv( Element& y, const Element& x ) const
		{
			Element one(0, 1);
			return quo(y,one,x);
		}

		Element& invin( Element& y ) const
		{
			Element x = y;
			return inv(y, x);
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

		static inline integer getMaxModulus()
		{ return NTL_zz_p::getMaxModulus(); }
		/** Write a description of the field */
		// Oustide of class definition so write(ostream&,const Element&) from
		// UnparametricField still works.

		std::ostream& write( std::ostream& os ) const
		{
			return os << "Polynomial ring using NTL::zz_pEX";
		}
		std::ostream& write( std::ostream& os, const Element& x) const
		{	return Father_t::write(os, x); }

		/** Conversion to scalar types doesn't make sense and should not be
		 * used.  Use getCoeff or leadCoeff to get the scalar values of
		 * specific coefficients, and then convert them using coeffField()
		 * if needed.
		 */
		template< class ANY >
		ANY& convert( ANY& x, const Element& y ) const
		{ return x; }

	protected:
		CoeffField _CField;
	}; // end of class NTL_zz_pEX




	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::zz_pEX> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::zz_pEX> {
	public:
		typedef NTL::zz_pEX Element;
		UnparametricRandIter<NTL::zz_pEX>(const NTL_zz_pEX & F ,
						 const size_t& size = 0,
						 const size_t& seed = 0
						) :
			_size(size), _seed(seed)
		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ(_seed));
		}

		UnparametricRandIter<NTL::zz_pEX>(const UnparametricRandIter<NTL::zz_pEX>& R) :
			_size(R._size), _seed(R._seed)

		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ(_seed));
		}

		Element& random (Element& x) const
		{
			NTL::random(x, 1);
			return x;
		}

	protected:
		size_t _size;
		size_t _seed;
	}; // class UnparametricRandIters

} // end of namespace LinBox

#endif // __LINBOX_field_ntl_lzz_pEX_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

