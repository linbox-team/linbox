/* linbox/field/ntl.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 * Copyright (C) 2011 LinBox
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file field/ntl/ntl-ZZ_p.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_zz_p_H
#define __LINBOX_field_ntl_zz_p_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <sys/time.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include <givaro/zring.h>
#include "linbox/field/field-traits.h"

#include "linbox/vector/blas-vector.h"


#include "linbox/integer.h"
#include "ntl-zz.h"

namespace Givaro
{


	/** Conversion of field element to an Integer.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * This done by converting to a std::string : inefficient but correct.
	 * @return reference to Integer.
	 * @param x reference to Integer to contain output (reference returned).
	 * @param y constant reference to field element.
	 */
	template <>
	Integer& Caster(Integer& x, const NTL::ZZ_p& y)
	{
		NTL::ZZ iy = y._ZZ_p__rep;

		long nb = NTL::NumBytes(iy);
		unsigned char *txt;
		typedef unsigned char u_char;
		txt = new u_char[nb + 68];
		//			   if (!txt) Error("out of memory");
		BytesFromZZ(txt, iy, nb);

		x = 0;
		for (ptrdiff_t i = 0; i < nb; i++) {
			x += Integer( txt[i] )<< int32_t(8*i) ;
		}
		delete [] txt;
		return x;
	}

	//dpritcha
	template<>
	double& Caster(double& x, const NTL::ZZ_p& y)
	{
		x = NTL::to_double(NTL::rep(y));
		return x;
	}

	/**\brief Initialization of field element from an Integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * This done by converting to a std::string : inefficient but correct.
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y Integer.
	 \ingroup field
	 */
	template <>
	NTL::ZZ_p& Caster(NTL::ZZ_p& x, const Integer& y)
	{
		return x = NTL::to_ZZ_p( NTL::to_ZZ( (static_cast<const std::string>(y)).c_str() ) );
	}
	template <>
	NTL::ZZ_p& Caster(NTL::ZZ_p& x, const double& y)
	{
		return x = NTL::to_ZZ_p( NTL::to_ZZ((long)(y) ) );
	}
	template <>
	NTL::ZZ_p& Caster(NTL::ZZ_p& x, const int& y)
	{
		return x = NTL::to_ZZ_p( NTL::to_ZZ((long)(y) ) );
	}
	template <>
	NTL::ZZ_p& Caster(NTL::ZZ_p& x, const unsigned long& y)
	{
		return x = NTL::to_ZZ_p( NTL::to_ZZ((long)(y) ) );
	}
	template <>
	NTL::ZZ_p& Caster(NTL::ZZ_p& x, const unsigned int& y)
	{
		return x = NTL::to_ZZ_p( NTL::to_ZZ((long)(y) ) );
	}
} // namespace Givaro


// Namespace in which all LinBox library code resides
namespace LinBox
{


	class NTL_ZZ_p_Initialiser {
	public :
		NTL_ZZ_p_Initialiser( const Integer & q, size_t e = 1) {
			linbox_check(e == 1);
			if ( q > 0 )
				NTL::ZZ_p::init(NTL::to_ZZ((std::string(q)).data())); // it's an error if q not prime, e not 1
		}

		template <class ElementInt>
		NTL_ZZ_p_Initialiser(const ElementInt& d) {
			NTL::ZZ_p::init (NTL::to_ZZ(d));
		}

		NTL_ZZ_p_Initialiser (const NTL::ZZ& d) {
			NTL::ZZ_p::init(d);
		}

		NTL_ZZ_p_Initialiser () { }

	};

            // CP: to be changed in to a Givaro::GeneralRingRandIter ?
        template <class Element>
	class UnparametricRandIter;
        /**
	 *
	 * @brief Wrapper of zz_p from NTL.
	 * Uses nice mod p via floating pt trick.
	 *
	 */
	struct NTL_ZZ_p: public NTL_ZZ_p_Initialiser, public Givaro::UnparametricOperations<NTL::ZZ_p> {
		typedef NTL::ZZ_p Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;

		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mOne ;

		/** @name NTL_ZZ_p
		 * @brief Arbitrary precision integers modulus a positive integer.

		 * While NTL allows any integer to serve as the modulus, only prime
		 * moduli yield fields.  Therefore, while arthmetic operations may be
		 * valid for any modulus, only prime moduli are supported in this
		 * implementation.  The primality of the modulus will not be checked, so
		 * it is the programmer's responsibility to supply a prime modulus.
		 * These specializations allow the \ref Givaro::ZRing template class to be
		 * used to wrap NTL's <code>ZZ_p</code> class as a LinBox field.
		 */
		//@{
		//! @param q,e
		NTL_ZZ_p(integer q, size_t e = 1) :
			NTL_ZZ_p_Initialiser(q,e),Father_t ()
			,zero( NTL::to_ZZ_p(0)),one( NTL::to_ZZ_p(1)),mOne(-one)
		{
			// no default - allow initialization of ZZ_p directly by user.
		}

		//! @param d,e
		NTL_ZZ_p( NTL::ZZ d, size_t e = 1) :
			NTL_ZZ_p_Initialiser(d),Father_t()
			,zero( NTL::to_ZZ_p(0)),one( NTL::to_ZZ_p(1)),mOne(-one)
		{
			linbox_check(e == 1);
		}

		//! NULL constructor
		NTL_ZZ_p() :
			NTL_ZZ_p_Initialiser(), Father_t()
			,zero( NTL::to_ZZ_p(0)),one( NTL::to_ZZ_p(1)),mOne(-one)
		{}
		//@}

		Element& init(Element& x, const integer& y) const
		{
			return Caster(x,y);
		}

		Element& init(Element& x, const double& y) const
		{
			double z = fmod(y,NTL::to_double(Element::modulus()));
			if (z > 0) z += 0.5;
			else z -= 0.5;
			return x = NTL::to_ZZ_p(static_cast<long>(z)); //rounds towards 0
		}

		/** Specialization for NTL::ZZ
		 *
		 * @return reference to field element.
		 * @param x field element to contain output (reference returned)
		 * @param y NTL::ZZ.
		 */
		Element& init(Element& x, const NTL::ZZ& y) const
		{
			return x = NTL::to_ZZ_p( y );
		}

		Element& init(Element& x) const
		{
			return x = NTL::to_ZZ_p( 0L );
		}

		Element& init(Element& x, const Element& y) const
		{
			return x = y ;
		}

		template <class ANY> //dpritcha--FIX
		Element& init(Element& x, const ANY& y) const
		{
			return x = NTL::to_ZZ_p((long)(y));
		}

		/** Specialization for NTL::ZZ.
		 *
		 * @return reference to  NTL::ZZ
		 * @param x  NTL::ZZ to contain output (reference returned).
		 * @param y constant reference to field element.
		 */
		NTL::ZZ& convert(NTL::ZZ& x, const Element& y) const
		{
			return x = y._ZZ_p__rep;
		}

		/** Conversion of field element to an integer.
		 * This function assumes the output field element x has already been
		 * constructed, but that it is not already initialized.
		 * This done by converting to a std::string : inefficient but correct.
		 * @return reference to integer.
		 * @param x reference to integer to contain output (reference returned).
		 * @param y constant reference to field element.
		 */
		integer& convert(integer& x, const Element& y) const
		{
			NTL::ZZ iy = y._ZZ_p__rep;

			long nb = NTL::NumBytes(iy);
			unsigned char *txt;
			typedef unsigned char u_char;
			txt = new u_char[nb + 68];
			//			   if (!txt) Error("out of memory");
			BytesFromZZ(txt, iy, nb);

			x = 0;
			for (ptrdiff_t i = 0; i < nb; i++) {
				x += LinBox::integer( txt[i] )<<int32_t(8*i) ;
			}
			delete [] txt;
			return x;
		};

		double& convert(double& x, const Element& y) const
		{
			x = NTL::to_double(NTL::rep(y));
			return x;
		}

		template <class ANY>
		ANY& convert(ANY& x, const Element& y) const
		{
			return x = (ANY)(rep(y));
		}

		static inline integer maxCardinality()
		{
			return integer( -1 );
		}

		Element& pow( Element& res, const Element& x, long exp ) const
		{
			NTL::power( res, x, exp );
			return res;
		}

		Element& powin( Element& x, long exp ) const
		{
			return x = NTL::power(x,exp);
		}

		/** Cardinality.
		 * Return integer representing cardinality of the field.
		 * Returns the modulus of the field, which should be prime.
		 * @return integer representing cardinality of the field
		 */
		template <typename Residu_t> Residu_t &cardinality(Residu_t& c) const
		{
			return characteristic(c);
		}

		integer cardinality() const
		{
			return characteristic();
		}

		/** Characteristic.
		 * Return integer representing characteristic of the field.
		 * Returns the modulus of the field, which should be prime.
		 * @return integer representing characteristic of the field.
		 */
		integer& characteristic(integer& c) const
		//FIXME we shouldn't go thru long here as p may be larger than that.
		// check if NTL has cast ZZp to gmp integers.
		{
			std::stringstream s;
			s << Element::modulus();
			s >> c;
			return c;
			//return c = static_cast<integer>(to_long(Element::modulus()));
		}

		size_t& characteristic(size_t & c) const
		{
			return c = (int64_t)to_long(Element::modulus());
		}

		integer characteristic() const
		{
			integer c;
			return characteristic(c);
		}

		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element&
		inv(Element& x, const Element& y) const
		{
			return x = NTL::inv(y);
		}

		/** Zero equality.
		 * Test if field element is equal to zero.
		 * This function assumes the field element has already been
		 * constructed and initialized.
		 * In this specialization, NTL's IsZero function is called.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field element.
		 */
		bool isZero(const Element& x) const
		{
			return static_cast<bool>(IsZero(x));
		}

		/** One equality.
		 * Test if field element is equal to one.
		 * This function assumes the field element has already been
		 * constructed and initialized.
		 * In this specialization, NTL's IsOne function is called.
		 * @return boolean true if equals one, false if not.
		 * @param  x field element.
		 */
		bool isOne(const Element& x) const
		{
			return static_cast<bool>(IsOne(x));
		}

		/** MOne equality.
		 * Test if field element is equal to one.
		 * This function assumes the field element has already been
		 * constructed and initialized.
		 * In this specialization, NTL's IsMOne function is called.
		 * @return boolean true if equals one, false if not.
		 * @param  x field element.
		 */
		bool isMOne(const Element& x) const
		{
		Element y ; neg(y,x);
			return isOne(y);
		}

            /** Unit test.
             * Test if field element is invertible.
             * This function assumes the field element has already been
             * constructed and initialized.
             * In this specialization, NTL's InvModStatus function is called.
             * inline long InvModStatus(ZZ& x, const ZZ& a, const ZZ& n)
             * // if gcd(a,n) = 1, then ReturnValue = 0, x = a^{-1} mod n
             * // otherwise, ReturnValue = 1, x = gcd(a, n)
             * @return boolean true if invertible, false if not.
             * @param  x field element.
             */
		bool isUnit(const Element& x) const
            {
                NTL::ZZ d;
                return !NTL::InvModStatus(d,rep(x),NTL::ZZ_p::modulus());
            }

		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes both field elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 */
		Element& invin(Element& x) const
		{
			return x = NTL::inv(x);
		}

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream& write(std::ostream& os) const
		{
			return os << "unparameterized field Element with p = "
			<< Element::modulus();
		}

		/** Print field.
		 * @return output stream to which field is written.
		 * @param os  output stream to which field is written.
		 * @param x
		 */
		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return Givaro::UnparametricOperations<Element>::write(os,x);
		}
	};

	template <class Ring>
	struct ClassifyRing;

	template <>
	struct ClassifyRing<NTL_ZZ_p > {
		typedef RingCategories::ModularTag categoryTag;
	};

	/// Constructor for random field element generator
	template <>
	class UnparametricRandIter<NTL::ZZ_p> {
	protected:
		integer _size;
        uint64_t _seed;
        const NTL_ZZ_p & _ring;

	public:
        typedef NTL::ZZ_p Element;
        typedef integer Residu_t;

		UnparametricRandIter(const NTL_ZZ_p & F,
                             const uint64_t seed = 0,
                             const Residu_t& size = 0) :
                _size(size), _seed(seed), _ring(F)
		{
                    if (_seed == 0) _seed = static_cast<uint64_t>(std::time(nullptr));

			integer cardinality;
			F.cardinality(cardinality);
			if (_size > cardinality)
				_size = 0;

#ifdef TRACE
			std::cout << "created random generator with size " << _size
			<< " and seed " << _seed << std::endl;
#endif // TRACE

			// Seed random number generator
			NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(_seed)));
		}

        const NTL_ZZ_p& ring() const { return _ring; }


		// UnparametricRandIter<NTL::ZZ_p>(const NTL_ZZ_p& R) :
			// _size(R._size), _seed(R._seed)
		// {
			// if(_seed == 0)
				// NTL::SetSeed(NTL::to_ZZ(time(nullptr)));
			// else
				// NTL::SetSeed(NTL::to_ZZ( static_cast<long>(_seed)) );
		// }


		/// Random field element creator.
		NTL::ZZ_p& random(NTL::ZZ_p& x) const
		{
			if (_size == 0) {
				return x = NTL::random_ZZ_p();
			}
			else {
				return x = NTL::to_ZZ_p(NTL::RandomBnd(Caster<NTL::ZZ,Residu_t>(_size)));
			}
		}


	};


} // namespace LinBox

#endif // __LINBOX_field_ntl_zz_p_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
