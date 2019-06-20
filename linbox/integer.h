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
 * ring for use in integer matrix computation, see <givaro/zring.h> or ring/ntl/ntl-zz.h.
 */

#ifndef __LINBOX_integer_H
#define __LINBOX_integer_H

#include <givaro/givconfig.h>
#include "linbox/linbox-config.h"
#include <gmp++/gmp++.h>


namespace LinBox
{
	/*! Integers in LinBox.
	 * Integer representation from <a href=http://ljk.imag.fr/CASYS/LOGICIELS/givaro/>Givaro</a>.
	 * @ingroup integers
	 */
	typedef Givaro::Integer integer;
	typedef Givaro::Integer Integer;

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
	inline double naturallog(const Givaro::Integer& a) {
		return Givaro::naturallog(a);
	}
}


#include <givaro/givcaster.h>



namespace LinBox { /*  signedness of integers */
	/*! Positiveness of an integer.
	 * Essentially usefull in debug mode to avoid compiler warnings
	 * about comparison always true for some unsigned type.
	 * @param x integer
	 * @return \c true iff \c x>=0.
	 */
	//@{
    //


    template< class T>
    inline typename std::enable_if<!std::is_unsigned<T>::value, bool>::value
    isPositive (const T& t) {
        return t >= 0;
    }

    template< class T>
    inline typename std::enable_if<std::is_unsigned<T>::value, bool>::value
    isPositive (const T& t) {
        return true;
    }

    template<class T>
    inline typename std::enable_if<!std::is_unsigned<T>::value, bool>::value
    isNegative (const T& t) {
        return (!isPositive<T>(t));
    }

    template<class T>
    inline typename std::enable_if<std::is_unsigned<T>::value, bool>::value
    isNegative (const T& t) {
        return false;
    }

	template<typename IntType>
	inline bool isOdd (const IntType& value)
	{
		return value & IntType(1);
	}

	template<typename IntType>
	inline bool isEven (const IntType& p)
	{
		return !isOdd(p);
	}

}

namespace LinBox
{ /*  indexDomain : used only once in linbox/algorithms/rational-solver.inl */


	/** Class used for permuting indices.
	 * For example, create a vector <code>(0 1 2 ...)</code> over \c
	 * size_t, then apply a permutation to it using a \c BlasMatrixDomain to
	 * get the natural representation of the permutation.
	 * @bug does not belong here
	 */
	class indexDomain
	{
	public:
		typedef size_t Element;
		typedef Element* Element_ptr ;
		typedef const Element* ConstElement_ptr ;

	public:
		typedef indexDomain Father_t;
		indexDomain() {};
		
		size_t init(size_t& dst) const {
			return dst = static_cast<size_t>(0);
		}
		
		template <class ANY>
		size_t init(size_t& dst, const ANY& src) const {
			return dst = static_cast<size_t>(src);
		}
		
		template <class ANY>
		size_t assign(ANY& dst, const size_t& src) const {
			return dst = static_cast<ANY>(src);
		}

		int characteristic() const { return 0 ; }
	};
}


#endif // __LINBOX_integer_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
