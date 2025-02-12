/* linbox/field/ntl-lzz_p.h
 * Copyright (C) 1999-2005 W. J. Turner,
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
 *
 */

/*! @file field/ntl/ntl-lzz_p.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_lzz_p_H
#define __LINBOX_field_ntl_lzz_p_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <time.h>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include <NTL/lzz_p.h>
#include <NTL/ZZ.h>

#include <givaro/zring.h>
#include <givaro/unparametric-operations.h>
#include "linbox/field/field-traits.h"
#include "linbox/integer.h"

namespace Givaro
{

        /** Initialization of field element from an integer.
         * This Uses NTL's <code>to_zz_p</code> function.
         *
         * @return reference to field element.
         * @param x field element to contain output (reference returned).
         * @param y Integer.
         */
	template <>
	NTL::zz_p& Caster(NTL::zz_p& x, const Integer& y)
	{
		return x = NTL::to_zz_p( static_cast<int64_t>(y)%NTL::zz_p::modulus() );
	}
	template <>
	NTL::zz_p& Caster(NTL::zz_p& x, const double& y)
	{
		return x = NTL::to_zz_p(long(y)%NTL::zz_p::modulus());
	}

        /** Conversion of field element to an Integer.
         * This function assumes the output field element x has already been
         * constructed, but that it is not already initialized.
         * For now, this is done by converting the element type to a C++
         * long and then to the Integer type through the use of static cast and
         * NTL's to_long function.
         * This, of course, assumes such static casts are possible.
         * This function should be changed in the future to avoid using long.
         * @return reference to Integer.
         * @param x reference to Integer to contain output (reference returned).
         * @param y constant reference to field element.
         */
    template <>
	Integer& Caster(Integer& x, const NTL::zz_p& y)
	{
		return x = int64_t(rep(y));
	}
} // namespace Givaro

// Namespace in which all LinBox library code resides
namespace LinBox
{

	class NTL_zz_p_Initialiser {
	public :
		NTL_zz_p_Initialiser( const Integer & q, size_t e = 1) {
			linbox_check(e == 1);
			if ( q > 0 )
                    // NTL::zz_p::init(NTL::to_ZZ((std::string(q)).data())); // it's an error if q not prime, e not 1
				NTL::zz_p::init(int64_t(q)); // it's an error if q not prime, e not 1
		}

		template <class ElementInt>
		NTL_zz_p_Initialiser(const ElementInt& d) {
			NTL::zz_p::init (NTL::to_ZZ(d));
		}

		NTL_zz_p_Initialiser () { }

	};

        /**
         * \brief long ints modulo a positive integer.
         *
         * While NTL allows any int to serve as the modulus, only prime
         * moduli yield fields.  The primality of the modulus will not be checked, so
         * it is the programmer's responsibility to supply a prime modulus if a field is
         * wanted.
         * These specializations allow the \ref Givaro::ZRing template class to be
         * used to wrap NTL's <code>zz_p</code> class as a LinBox field.
         * Uses nice trick for mod p via floating point.
         \ingroup field
        */
	struct NTL_zz_p: public NTL_zz_p_Initialiser, public Givaro::UnparametricOperations<NTL::zz_p> {
		typedef NTL::zz_p Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<NTL::zz_p> RandIter;

		const Element zero,one,mOne ;


            //public Givaro::ZRing<Element> {
		NTL_zz_p(integer p, size_t e = 1) :
                NTL_zz_p_Initialiser(p,e),Father_t ()
			,zero( NTL::to_zz_p(0)),one( NTL::to_zz_p(1)),mOne(-one)
            {}

		NTL_zz_p() :
                NTL_zz_p_Initialiser(2,1), Father_t()
			,zero( NTL::to_zz_p(0)),one( NTL::to_zz_p(1)),mOne(-one)
            {}

		Element &init (Element &x) const
            {
                return x = NTL::to_zz_p(0);
            }

		Element& init(Element& x, const double& y) const
            {
                double z = fmod(y,(double)Element::modulus());
                if (z > 0) z += 0.5;
                else z -= 0.5;
                return x = NTL::to_zz_p(static_cast<long>(z)); //rounds towards 0
            }

		Element &init (Element &x, const integer &y) const
            {
                NTL::ZZ tmp= NTL::to_ZZ(std::string(y).data());
                return x = NTL::to_zz_p(tmp);
            }

/*
  template <class ANY>
  Element& init(Element& x, const ANY& y) const
  {
  return x = NTL::to_zz_p((long)(y));
  }
*/

		Element& init(Element& x, const NTL::zz_p & y) const
            {
                return x = y ;
            }

		template <class ANY>
		ANY& convert(ANY& x, const Element& y) const
            {
                return Givaro::Caster(x,y);
            }

		static inline integer maxCardinality()
            {
                return integer( int64_t(NTL_SP_BOUND) );
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

		template <typename Residu_t> Residu_t &characteristic(Residu_t& c) const
            {
                return c = static_cast<int64_t>(Element::modulus());
            }

		integer characteristic() const
            {
                return static_cast<int64_t>(Element::modulus());
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
                return static_cast<bool>(NTL::IsZero(x));
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
                return static_cast<bool>(NTL::IsOne(x));
            }

            /** Unit test.
             * Test if field element is invertible.
             * This function assumes the field element has already been
             * constructed and initialized.
             * @return boolean true if invertible, false if not.
             * @param  x field element.
             */
		bool isUnit(const Element& x) const
            {
                long d,u;
                Givaro::extended_euclid(u,d,rep(x),NTL::zz_p::modulus());
                return (d==1) || (d==-1);
//                 return !NTL::InvModStatus(d,rep(x),NTL::zz_p::modulus());
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
                Element y = x ; negin(y);
                return isOne(y);
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
		std::ostream& write(std::ostream& os, const std::string& ) const
            {
                return os << "unparameterized field Element with p = "
                          << Element::modulus();
            }

		std::ostream &write (std::ostream &os, const Element &x) const { return Givaro::UnparametricOperations<Element>::write(os,x); }
	};

	template <class Ring>
	struct ClassifyRing;

	template <>
	struct ClassifyRing<NTL_zz_p> {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::zz_p> {
	public:
		typedef NTL::zz_p Element;
        typedef size_t Residu_t;

            /// Constructor for random field element generator

		UnparametricRandIter(const NTL_zz_p & F,
                             const uint64_t seed=0,
                             const Residu_t& size=0) :
                _size(size), _seed(seed), _ring(F)
            {
                if (_seed == 0) _seed = static_cast<uint64_t>(std::time(nullptr));

                size_t cardinality;
                F.cardinality(cardinality);
                if (_size > cardinality)
                    _size = 0;

#ifdef TRACE
                std::cout << "created random generator with size " << _size
                          << " and seed " << _seed << std::endl;
#endif // TRACE

                    // Seed random number generator
                std::stringstream s;
                NTL::ZZ x;
                s << _seed;
                s >> x;
                NTL::SetSeed(x);//NTL::to_ZZ(static_cast<long>(_seed)));
            }

            /// Random field element creator.
		NTL::zz_p& random(NTL::zz_p& x) const
            {
                if (_size == 0) {
                    return x = NTL::random_zz_p();
                }
                else {
                    return x = NTL::to_zz_p(NTL::RandomBnd(static_cast<int64_t>(_size)));
                }
            }
        const NTL_zz_p& ring() const { return _ring; }
	protected :
		size_t _size;
        uint64_t _seed ;
        const NTL_zz_p& _ring;
	};

} // namespace LinBox

#endif // __LINBOX_field_ntl_lzz_p_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
