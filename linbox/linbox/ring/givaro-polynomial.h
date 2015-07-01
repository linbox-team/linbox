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
#include "linbox/element/givaro-polynomial.h"


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

	template<class PolyCont>
	PolyCont& factor (PolyCont& factors, 
			  std::vector<unsigned long>& exp,
			  const Polynomial& P);
	
};

	
#ifdef __LINBOX_HAVE_NTL
}
#include "linbox/field/ntl-ZZ.h"
#include "NTL/ZZXFactoring.h"
namespace LinBox{
template <>
template <>
std::vector<GivPolynomial<integer>* >& 
GivPolynomialRing<UnparametricField<integer>,Dense>::factor (std::vector<GivPolynomial<integer>* >& factors, 
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
			factors[i] = new GivPolynomial<integer>( deg(ntlfactors[i].a)+1 );
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff(t,ntlfactors[i].a,j);
				NTLIntDom.convert( factors[i]->operator[](j), t );
			}
			exp[i] = ntlfactors[i].b;
		}
		return factors;
}

#include <linbox/field/PID-integer.h>
template <>
template <>
std::vector<GivPolynomial<integer>* >& 
GivPolynomialRing<PID_integer,Dense>::factor<std::vector<GivPolynomial<integer>* > > (std::vector<GivPolynomial<integer>* >& factors, 
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
			factors[i] = new GivPolynomial<integer>( deg(ntlfactors[i].a)+1 );
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff(t,ntlfactors[i].a,j);
				NTLIntDom.convert( factors[i]->operator[](j), t );
			}
			exp[i] = ntlfactors[i].b;
		}
		return factors;
}

template <>
template <>
std::vector<GivPolynomial<NTL::ZZ>* >& 
GivPolynomialRing< NTL_ZZ , Dense>::factor<std::vector<GivPolynomial<NTL::ZZ>* > > (std::vector<GivPolynomial<NTL::ZZ>* >& factors, 
										    std::vector<unsigned long>& exp,
										    const GivPolynomial<NTL::ZZ> &P)
{
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < P.size(); ++i){
			NTL::SetCoeff (f, i, P[i]);
		}
		NTL::vec_pair_ZZX_long ntlfactors;
		NTL::ZZ c;
		NTL::factor (c, ntlfactors, f);
			
		NTL::ZZ t; 
		NTL_ZZ NTLIntDom;
		factors.resize(ntlfactors.length());
		exp.resize(ntlfactors.length());
		for (int i= 0; i<ntlfactors.length(); ++i) {
			factors[i] = new GivPolynomial<NTL::ZZ>( deg(ntlfactors[i].a)+1 );
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff( factors[i]->operator[](j),ntlfactors[i].a,j);
			}
			exp[i] = ntlfactors[i].b;
		}
		return factors;
} 

#endif

template <>
template <>
std::vector<GivPolynomial<double> *>& 
GivPolynomialRing<Modular<double>,Dense>::factor (std::vector<GivPolynomial<double>* > & factors, 
						  std::vector<unsigned long>& exp,
						  const GivPolynomial<double>& P)
{
	integer charac;
	_domain.characteristic(charac);
	double p = charac;
	Poly1FactorDom<Modular<double>,Dense> PFD(*this);
	std::vector<givvector<double> > factors2;
	PFD.CZfactor ( factors2, exp, static_cast<givvector<double> >(P),p);

	//std::cerr<<"factorization done"<<std::endl;
	factors.resize(factors2.size());
	std::vector<GivPolynomial<double>* >::iterator itf = factors.begin();
	std::vector<givvector<double> >::const_iterator itf2 = factors2.begin();
	for (; itf2 != factors2.end();++itf,++itf2){
		*itf = new GivPolynomial<double>(*itf2);
		//std::cerr<<"converting factor"<<(*itf)<<std::endl;
		for (size_t i=0; i< (*itf)->size();++i)
			_domain.divin((*itf)->operator[](i),(*itf)->operator[]((*itf)->size()-1));
		_domain.assign((*itf)->operator[]((*itf)->size()-1),1.0);
	}
	return factors;
}	

} // namespace LinBox


#endif // __GIVAROPOLYNOMIAL_H
