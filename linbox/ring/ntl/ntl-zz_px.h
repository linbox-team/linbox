/* Copyright (C) 2005 LinBox
 * written by Daniel Roche, August 2005
 * Copyright (C) 2011 LinBox
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

/*! @file ring/ntl/ntl-ZZ_p.h
 * @ingroup ring
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_zz_px_H
#define __LINBOX_field_ntl_zz_px_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <vector>
#include <NTL/ZZ_pX.h>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include <givaro/zring.h>
#include "linbox/ring/ntl/ntl-zz_p.h"
#include "linbox/integer.h"


// namespace in which all LinBox code resides
namespace LinBox
{
	class NTL_ZZ_pX_Initialiser {
	public :
		NTL_ZZ_pX_Initialiser( const Integer & q, size_t e = 1) {
			// linbox_check(e == 1);
			if ( q > 0 )
				NTL::ZZ_p::init(NTL::to_ZZ((std::string(q)).data())); // it's an error if q not prime, e not 1
		}

		// template <class ElementInt>
		// NTL_ZZ_p_Initialiser(const ElementInt& d) {
		// NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		// NTL_ZZ_p_Initialiser (const NTL::ZZ& d) {
		// NTL::ZZ_p::init(d);
		// }

	};

	/** Ring (in fact, a unique factorization domain) of polynomial with
	 * coefficients in class NTL_ZZ_p (integers mod a wordsize prime).
	 * All the same functions as any other ring, with the addition of:
	 * Coeff (type), CoeffField (type), getCoeffField, setCoeff, getCoeff,
	 * leadCoeff, deg
	 */
	class NTL_ZZ_pX : public NTL_ZZ_pX_Initialiser, public Givaro::UnparametricOperations<NTL::ZZ_pX> {
	public:
		typedef NTL::ZZ_pX Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mOne ;


		typedef NTL_ZZ_p CoeffField;
		typedef NTL::ZZ_p Coeff;
		// typedef NTL::ZZ_pX Element;

		/** Standard LinBox field constructor.  The paramters here
		 * (prime, exponent) are only used to initialize the coefficient field.
		 */
		NTL_ZZ_pX( const integer& p, size_t e = 1 ) :
			NTL_ZZ_pX_Initialiser(p,e),Father_t ()
			, zero( NTL::to_ZZ_pX(0)),one( NTL::to_ZZ_pX(1)),mOne(-one)
			, _CField(p,e)
		{}

		/** Constructor from a coefficient field */
		NTL_ZZ_pX( CoeffField cf ) :
			NTL_ZZ_pX_Initialiser(cf.cardinality()),Father_t ()
			,zero( NTL::to_ZZ_pX(0)),one( NTL::to_ZZ_pX(1)),mOne(-one)
			,_CField(cf)
		{}

		/** Initialize p to 0 */
		Element& init( Element& p ) const
		{
			return p = 0;
		}

		/** Initialize p to the constant y (p = y*x^0) */
		template <class ANY>
		Element& init( Element& p, const ANY& y ) const
		{
			Coeff temp;
			_CField.init( temp, y );
			return p = temp;
		}

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

		bool isUnit(const Element& x) const
		{
			return ( (this->deg(x) == 0) &&
				 ( _CField.isUnit( NTL::ConstTerm(x) ) ) );
		}

		bool isMOne (const Element& x) const
		{
			return ( (this->deg(x) == 0) &&
				 ( _CField.isMOne( NTL::ConstTerm(x) ) ) );
		}

		/** The LinBox field for coefficients */
		const CoeffField& getCoeffField() const
		{
			return _CField;
		}

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

		/** Get characteristic of the field - same as characteristic of
		 * coefficient field. */
		integer& characteristic( integer& c ) const
		{
			return _CField.characteristic(c);
		}

		/** Get the cardinality of the field.  Since the cardinality is
		 * infinite, by convention we return -1.
		 */
		integer& cardinality( integer& c ) const
		{
			return c = static_cast<integer>(-1);
		}

		static inline integer maxCardinality()
		{
			return CoeffField::maxCardinality();
		}

		/** Write a description of the field */
		// Oustide of class definition so write(ostream&,const Element&) from
		// Givaro::ZRing still works.
		std::ostream& write( std::ostream& os ) const
		{
			return os << "Polynomial ring using NTL::ZZ_pX";
		}

		std::ostream &write (std::ostream &os, const Element &x) const { return Givaro::UnparametricOperations<Element>::write(os,x); }


	private:
		/** Conversion to scalar types doesn't make sense and should not be
		 * used.  Use getCoeff or leadCoeff to get the scalar values of
		 * specific coefficients, and then convert them using coeffField()
		 * if needed.
		 */
		template< class ANY >
		ANY& convert( ANY& x, const Element& y ) const
		{
			return x;
		}

		CoeffField _CField;
	}; // end of class NTL_ZZ_pX


} // end of namespace LinBox

#endif // __LINBOX_field_ntl_zz_px_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
