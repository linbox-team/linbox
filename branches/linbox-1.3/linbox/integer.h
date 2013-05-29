/* Copyright(c)'94-97 by Givaro Team
 * Copyright(c)'2000-2002 by LinBox Team
 *  ========LICENCE========
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
 * Created by M. Samama, T. Gautier
 *
 * Modified Jean-Guillaume.Dumas <Jean-Guillaume.Dumas@imag.fr>
 *          B. David Saunders <saunders@cis.udel.edu>,
 *          Bradford Hovinen <hovinen@cis.udel.edu>
 *          Gilles Villard <Gilles.Villard@ens-lyon.fr>
 *                        JGD Random functions back.
 *                        (2002/02/12 16:05:24)
 *
 */

/** @file integer.h
 * \ingroup linbox
 * \brief This is a representation of arbitrary integers.
 *
 * It is a wrapper of <a href=http://gmplib.org>GMP</a> integers.  Arithmetic operations are via
 * \c C++ infix operator forms (eg. \c a*b) . It is for ``casual'' uses such as characteristics and
 * cardinalities and when initializing field elements.  The integers are also represented as a
 * LinBox ring for use in integer matrix computation, see PID-integer.h  or see  field/ntl-ZZ.h.
 */

#ifndef __LINBOX_integer_H
#define __LINBOX_integer_H

#ifndef INT32_MAX
#define INT32_MAX (2147483647L)
#endif

//#include <cstdint>
#include "linbox/linbox-config.h"
#include "givaro/givconfig.h"
#include "gmp++/gmp++.h"
#include <cfloat> // BB : needed on some rare platforms...

using std::ptrdiff_t;


#ifdef __LINBOX_HAVE_STDINT_H
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
// else ??
#endif
#include <stdint.h>
#ifndef INT32_MAX
#error "INT32_MAX is not defined. It should at least be defined in Givaro..."
#endif
#endif


#ifndef FFLAFLAS_VERSION
#define FFLAFLAS_VERSION __LINBOX_FFLAFFLAS_VERSION
#endif


namespace LinBox
{
	/*! Integers in LinBox.
	 * Integer representation from <a href=http://ljk.imag.fr/CASYS/LOGICIELS/givaro/>Givaro</a>.
	 * @ingroup integers
	 */
	typedef Givaro::Integer integer;
	typedef Givaro::Integer Integer;

#if 0
// These integer types are defined by Givaro
	/// int8_t.
	typedef signed __LINBOX_INT8 int8_t;
	/// int16_t.
	typedef signed __LINBOX_INT16 int16_t;

	/** @brief This is a representation of 32 bit ints, usually equivalent to \c int.
	 *
	 * The use of \c int32_t ensures you are working with
	 * 32 bit signed ints, \f$[-2^{31}\dots2^{31})\f$.  Similarly, \ref int8_t, \ref int16_t, and \ref int64_t are defined.
	 */
	typedef signed __LINBOX_INT32 int32_t;

	/// int64_t.
	typedef signed __LINBOX_INT64 int64_t;

	/// unsigned int8_t.
	typedef unsigned __LINBOX_INT8 uint8_t;
	/// unsigned int16_t.
	typedef unsigned __LINBOX_INT16 uint16_t;

	/** This is a representation of 32 bit unsigned ints, usually
	 * equivalent to `<code>unsigned int</code>'.
	 *
	 * The use of `uint32_t' ensures you are working with
	 * 32 bit unsigned ints, \f$[0\cdots 2^32[\f$.  Similarly, uint8_t, uint16_t, and uint64_t are defined.
	 */
	typedef unsigned __LINBOX_INT32 uint32_t;

	/// unsigned int64_t.
	typedef unsigned __LINBOX_INT64 uint64_t;

#endif

	// Huh? -bds
	template< class T >
	T abs( const T& a ) { return( a <= 0 ? a * -1 : a ); }

} // LinBox namespace


// Dependency to GIVARO >= 3.7.2
#include <givaro/givspyinteger.h>
namespace LinBox
{

    /*! @internal
	 * Spy structure to have access to protected members of Givaro::Integer.
	 */
    using Givaro::SpyInteger;

} // LinBox namespace




// Dependency to GIVARO >= 3.3.4
/* givaro/givconfig.h so provides the fixed width integer types such as
 * int16_t, uint8_t, etc.  The typenames int16, uint8, etc are no longer used
 * in LinBox or Givaro.
 */
#include <givaro/givconfig.h>
#include <math.h>

#ifndef GIVARO_VERSION
#error "Givaro didn't tell us about his version !"
#endif


namespace LinBox
{

	/** Natural logarithm (ln).
	 * log(2) being close to 0.69314718055994531
	 * @param a integer.
	 * @return  ln(a).
	 */
#if (GIVARO_VERSION < 30305)
	inline double naturallog(const Givaro::Integer& a) {
		signed long int exp;
		double d = (double)mpz_get_d_2exp( &exp, (mpz_srcptr)(LinBox::SpyInteger::get_rep(a) ) );
		return (double)exp*0.69314718055994531+log(d);
	}
#else
	inline double naturallog(const Givaro::Integer& a) {
		return Givaro::naturallog(a);
	}
#endif
}

namespace LinBox { /*  signedness of integers */
	/*! Positiveness of an integer.
	 * Essentially usefull in debug mode to avoid compiler warnings
	 * about comparison always true for some unsigned type.
	 * @param x integer
	 * @return \c true iff \c x>=0.
	 */
	//@{
	template<class T>
	inline bool isPositive( const T & x) {
		return x>=0 ;
	}
	template<>
	inline bool isPositive(const uint8_t &) {
		return true ;
	}
	template<>
	inline bool isPositive(const uint16_t &) {
		return true ;
	}
	template<>
	inline bool isPositive(const uint32_t &) {
		return true ;
	}
#ifdef __APPLE__
	template<>
	inline bool isPositive(const unsigned long&) {
		return true ;
	}
#endif
	template<>
	inline bool isPositive(const uint64_t &) {
		return true ;
	}
	//@}
}

#if (GIVARO_VERSION < 30601)
namespace Givaro {
	template <typename Target, typename Source>
	Target& Caster (Target& t, const Source& s) {
		return t = static_cast<Target>(s);
	}
}
#else
#include <givaro/givcaster.h>
#endif

#ifdef GIVARO_USES_OPENMP // _OPENMP or others are present
#define LINBOX_USES_OPENMP 1
#endif

namespace LinBox {

	template<class U>
	inline bool IsNegative(const U & p)
	{
		return (p<0);
	}

	// or use integer_traits<T>::is_unsigned ??
	template<>
	inline bool IsNegative(const uint8_t & p)
	{
		return false;
	}

	template<>
	inline bool IsNegative(const uint16_t & p)
	{
		return false;
	}

	template<>
	inline bool IsNegative(const uint32_t & p)
	{
		return false;
	}

	template<>
	inline bool IsNegative(const uint64_t & p)
	{
		return false;
	}

}

#endif // __LINBOX_integer_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
