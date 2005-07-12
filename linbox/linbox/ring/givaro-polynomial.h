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
#include "givaro/givpoly1factor.h"
#include "linbox/integer.h"
#include "linbox/field/unparametric.h"
#include "linbox/field/ntl-ZZ.h"
#include "linbox/element/givaro-polynomial.h"
#include "NTL/ZZXFactoring.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	
	/** \brief polynomials with coefficients modulo some power of two
	\ingroup ring
	 *
	 * @param Polynomial type, e.g. std::vector<Field::Element>
	 */
template <class Domain, class StorageTag>
class GivPolynomialRing : public Poly1Dom<Domain,StorageTag>
{
public:

	//	using Poly1Dom<Domain,StorageTag>::eval;
	typedef GivPolynomial<typename Domain::Element> Element;

	typedef Element Polynomial;

	GivPolynomialRing () {}

	GivPolynomialRing (const Domain& D)
		: Poly1Dom<Domain,StorageTag>(D, Indeter()){}

	template<template< class >class Container>
	Container<Polynomial>& factor (Container<Polynomial>& factors, 
				       std::vector<unsigned long>& exp,
				       const Polynomial& P);
	
};

	
	//template<template< class >class Container>
template <>
std::vector<GivPolynomial<integer> >& 
GivPolynomialRing<UnparametricField<integer>,Dense>::factor (std::vector<GivPolynomial<integer> >& factors, 
							     std::vector<unsigned long>& exp,
							     const GivPolynomial<integer> &P)
{
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < P.size(); ++i){
			NTL::SetCoeff (f, i, NTL::to_ZZ((std::string( P[i] )).c_str()) );
		}
		NTL::vec_pair_ZZX_long ntlfactors;
		NTL::ZZ c;
		NTL::factor (c, ntlfactors, f);
			
		NTL::ZZ t; 
		NTL_ZZ NTLIntDom;
		factors.resize(ntlfactors.length());
		exp.resize(ntlfactors.length());
		for (int i= 0; i<ntlfactors.length(); ++i) {
			factors[i].resize( deg(ntlfactors[i].a)+1 );
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff(t,ntlfactors[i].a,j);
				NTLIntDom.convert( factors[i][j], t );
			}
			exp[i] = ntlfactors[i].b;
		}
		return factors;
}

template <>
std::vector<GivPolynomial<double> >& 
GivPolynomialRing<Modular<double>,Dense>::factor (std::vector<GivPolynomial<double> > & factors, 
						  std::vector<unsigned long>& exp,
						  const GivPolynomial<double>& P)
{
	integer charac;
	_domain.characteristic(charac);
	double p = charac;
	Poly1FactorDom<Modular<double>,Dense> PFD(*this);
	std::vector<givvector<double> > factors2;
	PFD.CZfactor ( factors2, exp, static_cast<givvector<double> >(P),p);

	factors.resize(factors2.size());
	std::vector<GivPolynomial<double> >::iterator itf = factors.begin();
	std::vector<givvector<double> >::const_iterator itf2 = factors2.begin();
	for (; itf2 != factors2.end();++itf,++itf2){
		*itf = *itf2;
		for (size_t i=0; i< itf->size();++i)
			_domain.divin((*itf)[i],(*itf)[itf->size()-1]);
		_domain.assign((*itf)[itf->size()-1],1.0);
	}
	return factors;
}	

} // namespace LinBox


#endif // __GIVAROPOLYNOMIAL_H
