/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/element/givaro-polynomial.h
 * Written by 
 * Clement Pernet
 *
 * See COPYING for license information.
 */

#ifndef __GIVAROPOLYNOMIAL_ELT_H
#define __GIVAROPOLYNOMIAL_ELT_H

#include <iostream>
#include "givaro/givpoly1.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	
	/** Polynomials over a domain
	 *
	 * @param Type of coefficients
	 */
template <typename T>
class GivPolynomial : public givvector<T>
{
public:
  
	GivPolynomial () : givvector<T>() {}
	
	GivPolynomial (size_t s) : givvector<T>(s) {}
	GivPolynomial (const std::vector<T>& p) : givvector<T>(p, givWithCopy()) {}

	template<typename X>
	struct rebind
	{
		typedef GivPolynomial<X> other;
	};
};

} // namespace LinBox


#endif // __GIVAROPOLYNOMIAL_ELT_H
