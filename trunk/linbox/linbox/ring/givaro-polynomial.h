/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/ring/givaro-polynomial.h
 * Written by 
 * Clement Pernet
 *
 * See COPYING for license information.
 */

#ifndef __GIVAROPOLYNOMIAL_H
#define __GIVAROPOLYNOMIAL_H

#include <iostream>
#include "givaro/givpoly1.h"
#include "linbox/integer.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	
	/** Ring of polynomials suing elements modulo some power of two
	 *
	 * @param Polynomial type, e.g. std::vector<Field::Element>
	 */
template <class Domain, class StorageTag>
class GivPolynomialRing : public Poly1Dom<Domain,StorageTag>
{
public:
	
	GivPolynomialRing (const Domain& d) : Poly1Dom(d){}

};
	

} // namespace LinBox


#endif // __GIVAROPOLYNOMIAL_H
