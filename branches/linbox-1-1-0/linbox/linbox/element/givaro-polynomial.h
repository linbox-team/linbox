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
	
	/** \brief Polynomials over a domain
	 *
	 * @param Type of coefficients
\ingroup element
	 */
template <typename T>
class GivPolynomial : public givvector<T>
{
	typedef GivPolynomial<T> Self_t;
public:
  
	GivPolynomial () : givvector<T>() {}
	
	GivPolynomial (size_t s) : givvector<T>(s) {}
	GivPolynomial (const std::vector<T>& p) : givvector<T>(p, givWithCopy()) {}
	GivPolynomial (const givvector<T>& p) : givvector<T>(p, givWithCopy()) {}

	template<typename X>
	struct rebind
	{
		typedef GivPolynomial<typename X::Element> other;
		
		void operator() (other *& P2, 
				 const Self_t& P1, 
				 const X& F)
		{
			P2 = new other(P1.size());
			typename Self_t::const_iterator it1 = P1.begin();
			typename other::iterator it2 = P2->begin();
			for (; it1 != P1.end(); ++it1, ++it2)
				F.init (*it2, *it1);
//  			for (size_t i=0; i < P1.size(); ++i)
//  				F.init ( (*P2)[i], P1[i]);
		}
	};
};

} // namespace LinBox


#endif // __GIVAROPOLYNOMIAL_ELT_H
