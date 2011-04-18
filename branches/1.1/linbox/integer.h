/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright(c)'94-97 by Givaro Team
 * Copyright(c)'2000-2002 by LinBox Team
 * see the COPYING file for license information.
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

//#include <cstdint>
#include "linbox/linbox-config.h"

#include "gmp++/gmp++.h"

#include <cfloat> // BB : needed on some rare platforms...
#ifdef __LINBOX_HAVE_STDINT_H
#define __STDC_LIMIT_MACROS
#include <stdint.h>
#endif

namespace LinBox
{
	/*! Integers in LinBox.
	 * Integer representation from <a href=http://ljk.imag.fr/CASYS/LOGICIELS/givaro/>Givaro</a>.
	 * @ingroup integers
	 */
	typedef Integer integer;

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



        /*! @internal
	 * Spy structure to have access to protected members of Givaro::Integer.
	 */
	struct SpyInteger
	{

	    struct InHeritsInteger : public integer {
	    protected:
	        friend struct SpyInteger;
	    };

	    static const InHeritsInteger::Rep* get_rep(const integer& i) {
        	return static_cast<const InHeritsInteger&>(i).get_rep();
	    }

	    static mpz_ptr get_mpz(integer& i) {
	        return static_cast<InHeritsInteger&>(i).get_mpz();
	    }
	    static mpz_ptr get_mpz(const integer& i) {
	        return const_cast<InHeritsInteger&>(static_cast<const InHeritsInteger&>(i)).get_mpz();
	    }
        };


}

// Dependency to GIVARO >= 3.3.4
/* givaro/givconfig.h so provides the fixed width integer types such as 
 * int16_t, uint8_t, etc.  The typenames int16, uint8, etc are no longer used 
 * in LinBox or Givaro.
 */
#include <givaro/givconfig.h>
#include <math.h>
// Natural logarithm of a
// log(2) being close to 0.69314718055994531
inline double naturallog(const Integer& a) {
  signed long int exp;
  double d = mpz_get_d_2exp( &exp, (mpz_ptr)(LinBox::SpyInteger::get_rep(a) ) );
  return (double)exp*0.69314718055994531+log(d);
}
#endif // __LINBOX_integer_H

