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
#include "NTL/ZZXFactoring.h"

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
	using Poly1Dom<Domain,StorageTag>:: Element;

	typedef typename Poly1Dom<Domain,StorageTag>::Element Polynomial;

	GivPolynomialRing () {}

	GivPolynomialRing (const Domain& D)
		: Poly1Dom<Domain,StorageTag>(D, Indeter()){}

	
	
	template<template< class >class Container>
	Container<Polynomial>& factor (Container<Polynomial>& factors, 
				       const Polynomial P) const
	{
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < P.size(); ++i)
			NTL::SetCoeff (f, i, NTL::to_ZZ((std::string( P[i] )).c_str()) );
		NTL::vec_pair_ZZX_long ntlfactors;
		NTL::ZZ c;
		NTL::factor (c, ntlfactors, f);
			
		NTL::ZZ t; 
		NTL_ZZ NTLIntDom;
		for (int i= 0; i<ntlfactors.length(); ++i) {
			factors[i].resize( deg(ntlfactors[i].a)+1 );
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff(t,ntlfactors[i].a,j);
				NTLIntDom.convert( factors[i][j], t );
			}
		}
		return factors;
	}

};
	

} // namespace LinBox


#endif // __GIVAROPOLYNOMIAL_H
