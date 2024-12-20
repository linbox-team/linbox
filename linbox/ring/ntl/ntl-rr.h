/* linbox/field/ntl-RR.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 * Copyright (C) 2011 LinBox
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
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

/*! @file field/ntl/ntl-RR.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_rr_H
#define __LINBOX_field_ntl_rr_H

#include <sstream>
#include "linbox/linbox-config.h"

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif


#include <NTL/tools.h>
#include <NTL/RR.h>

#include "linbox/util/debug.h"

#include <givaro/zring.h>
#include <givaro/unparametric-operations.h>
#include "linbox/field/field-traits.h"

#include "linbox/integer.h"

namespace Givaro
{
	/** Initialization of field element from an integer.
	 * Behaves like C++ allocator construct.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the integer type to a C++
	 * long and then to the element type through the use of static cast and
	 * NTL's to_RR function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to field element.
	 * @param x field element to contain output (reference returned).
	 * @param y integer.
	 */
	template <>
	NTL::RR& Caster(NTL::RR& x, const Integer& y)
	{
		std::stringstream s;
		s << y;
		s >> x;
		return x;
		//return x = NTL::to_RR(static_cast<long>(y));
	}
	template <> NTL::RR& Caster(NTL::RR& x, const double& y) { return x = NTL::to_RR((long)(y)); }

	template <> NTL::RR& Caster(NTL::RR& x, const int32_t& y) { return x = NTL::to_RR((long)(y)); }

	template <> NTL::RR& Caster(NTL::RR& x, const int64_t& y) { return x = NTL::to_RR((long)(y)); }

	template <> NTL::RR& Caster(NTL::RR& x, const uint32_t& y) { return x = NTL::to_RR((unsigned long)(y)); }

	template <> NTL::RR& Caster(NTL::RR& x, const uint64_t& y) { return x = NTL::to_RR((unsigned long)(y)); }

	/** Conversion of field element to an integer.
	 * This function assumes the output field element x has already been
	 * constructed, but that it is not already initialized.
	 * For now, this is done by converting the element type to a C++
	 * long and then to the integer type through the use of static cast and
	 * NTL's to_long function.
	 * This, of course, assumes such static casts are possible.
	 * This function should be changed in the future to avoid using long.
	 * @return reference to integer.
	 * @param x reference to integer to contain output (reference returned).
	 * @param y constant reference to field element.
	 */
	template <>
	Integer& Caster(Integer& x, const NTL::RR& y)
	{
		std::stringstream s;
		s << y;
		s >> x;
		return x;
		//return x = static_cast<Integer>(to_long(y));
	}
} // namespace Givaro




// Namespace in which all LinBox library code resides
namespace LinBox
{

	class NTL_RR_Initialiser {
	public :
		NTL_RR_Initialiser () { }
	};

    template<class XXX> class UnparametricRandIter;


	/** @name class RR.
	  \brief
	 * Rational number field.
	 * This field is provided as a convenience in a few places.
	 * Use with caution because expression swell.
	 *
	 * This specialization allows the \ref Givaro::ZRing template class to be
	 * used to wrap NTL's RR class as a LinBox field.
	 \ingroup field
	 */
	//@{


	struct NTL_RR: public NTL_RR_Initialiser, public Givaro::UnparametricOperations<NTL::RR> {
		typedef NTL::RR Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mOne ;

		NTL_RR() :
			NTL_RR_Initialiser(),Father_t ()
			,zero( NTL::to_RR(0)),one( NTL::to_RR(1)),mOne(-one)
		{
			// no default - allow initialization of ZZ_p directly by user.
		}


		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element& inv(Element& x, const Element& y) const
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

		/** invertibility.
		 * Test if field element is invertible.
		 * This function assumes the field element has already been
		 * constructed and initialized.
		 * @return boolean true if invertible, false if not.
		 * @param  x field element.
		 */
		bool isUnit(const Element& x) const
		{
			return !isZero(x);
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
			return os << "unparameterized field Element";
		}

		integer & cardinality(integer &c) const
		{
			return c = -1;
		}

		integer & characteristic(integer &c) const
		{
			return c = 0;
		}

		std::ostream &write (std::ostream &os, const Element &x) const {
			return Givaro::UnparametricOperations<Element>::write(os,x);
		}

		Element& init (Element& x) const
		{
			return x;
		}

		template <typename Src>
		Element& init (Element& x, const Src& s) const
		{
			return Caster (x, s);
		}

		template <typename T>
		T& convert (T &x, const Element &y) const
		{
			return Caster (x,y);
		}


	};

	template <>
	class UnparametricRandIter<NTL::RR> {

	protected:
		integer _size;
        uint64_t _seed;
        const NTL_RR& _ring;
	public:
		typedef NTL::RR Element ;
        typedef integer Residu_t;

		UnparametricRandIter(const NTL_RR & F,
                             const uint64_t seed = 0,
                             const Residu_t& size = 0) :
                _size(size), _seed(seed), _ring(F)
		{
                    if (_seed == 0) _seed = static_cast<uint64_t>(std::time(nullptr));

			// integer cardinality;
			// F.cardinality(cardinality);
			// if (_size > cardinality)
				// _size = 0;

#ifdef TRACE
			std::cout << "created random generator with size " << _size
			<< " and seed " << _seed << std::endl;
#endif // TRACE

			// Seed random number generator
			std::stringstream s; NTL::ZZ x;
			s << _seed;
			s >> x;
			NTL::SetSeed(NTL::to_ZZ(x));
		}

        const NTL_RR& ring() const { return _ring; }

		/** Random field element creator.
		 * This returns a random field element from the information supplied
		 * at the creation of the generator.
		 * This generator uses the built-in C++ random number generator instead of
		 * NTL's random function because the NTL function does not allow as much
		 * control over the sampling size as the generic LinBox template.  This
		 * specialization is included only to allow conversion to an NTL
		 * object.
		 * @return random field element
		 */
		Element& random(Element &elt) const
		{
			// NTL::random(elt);
			// Create new random elements
			if (_size == 0)
				elt = rand();
			else
				elt = static_cast<double>((double(rand())/RAND_MAX)*double(_size));

#ifdef TRACE
			double temp ;
			NTL::conv(temp, elt);
			std::cout << "random double = " << temp << "    random Element = " << elt << std::endl;
#endif // TRACE

			return elt;

		} // element& operator() (void)

	};

	template <class Ring>
	struct ClassifyRing;

	template<class Element>
	struct ClassifyRing<Givaro::ZRing<Element> >;

	template<>
	struct ClassifyRing<NTL_RR >{
		typedef RingCategories::ModularTag categoryTag;
	};
	//@}
} // namespace LinBox

#endif // __LINBOX_field_ntl_rr_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
