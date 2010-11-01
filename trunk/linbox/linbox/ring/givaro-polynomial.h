/* linbox/ring/givaro-polynomial.h
 * Copyright(C) LinBox
 * Written by 
 * Clement Pernet
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_givaropolynomial_H
#define __LINBOX_givaropolynomial_H

#include <iostream>
#include <givaro/givpoly1.h>
#include <givaro/givpoly1factor.h>
#include "linbox/integer.h"
#include "linbox/field/unparametric.h"
#include "linbox/field/givaro-field.h"
//#include "linbox/element/givaro-polynomial.h"


// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	
	/** \brief polynomials 
	\ingroup ring
	 *
	 * @param Polynomial type, e.g. std::vector<Field::Element>
	 */
template <class Domain, class StorageTag=Dense>
class GivPolynomialRing : public Poly1Dom<GivaroField<Domain>,StorageTag>
{
public:

	//	using Poly1Dom<Domain,StorageTag>::eval;
	typedef typename Poly1Dom<Domain,StorageTag>::Element Element;

	typedef Element Polynomial;

	GivPolynomialRing () {}

	GivPolynomialRing (const Domain& D)
		: Poly1Dom<GivaroField<Domain>,StorageTag>(D, Indeter()){}

	GivPolynomialRing (const Domain& D, const Indeter& I)
		: Poly1Dom<GivaroField<Domain>,StorageTag>(D, I){}

	template<class PolyCont>
	PolyCont& factor (PolyCont& factors, 
			  std::vector<unsigned long>& exp,
			  const Polynomial& P);
	
};

	
#ifdef __LINBOX_HAVE_NTL
}
#include "linbox/field/ntl-ZZ.h"
#include "NTL/ZZXFactoring.h"
namespace LinBox 
{
    typedef GivPolynomialRing<UnparametricField<integer>,Dense> GivPolIntDense;

    template <>
    template <>
    std::vector<GivPolIntDense::Element* >& 
    GivPolIntDense::factor (std::vector<GivPolIntDense::Element* >& factors, 
                            std::vector<unsigned long>& exp,
                            const GivPolIntDense::Element &P) {
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
            factors[i] = new GivPolIntDense::Element( deg(ntlfactors[i].a)+1 );
            for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
                NTL::GetCoeff(t,ntlfactors[i].a,j);
                NTLIntDom.convert( factors[i]->operator[](j), t );
            }
            exp[i] = ntlfactors[i].b;
        }
        return factors;
    }
}

#include <linbox/field/PID-integer.h>
namespace LinBox 
{
    typedef GivPolynomialRing<PID_integer,Dense> GivPolPIDIntDense;
    template <>
    template <>
    std::vector<GivPolPIDIntDense::Element* >& 
    GivPolPIDIntDense::factor<std::vector<GivPolPIDIntDense::Element* > > 
    (std::vector<GivPolPIDIntDense::Element* >& factors, 
     std::vector<unsigned long>& exp,
     const GivPolPIDIntDense::Element &P)
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
            factors[i] = new GivPolPIDIntDense::Element( deg(ntlfactors[i].a)+1 );
            for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
                NTL::GetCoeff(t,ntlfactors[i].a,j);
                NTLIntDom.convert( factors[i]->operator[](j), t );
            }
            exp[i] = ntlfactors[i].b;
        }
        return factors;
    }

    typedef GivPolynomialRing< NTL_ZZ , Dense> GivPolZZDense;
    template <>
    template <>
    std::vector<GivPolZZDense::Element* >& 
    GivPolZZDense::factor<std::vector<GivPolZZDense::Element* > > 
    (std::vector<GivPolZZDense::Element* >& factors, 
     std::vector<unsigned long>& exp,
     const GivPolZZDense::Element &P)
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
        factors.resize(ntlfactors.length());
        exp.resize(ntlfactors.length());
        for (int i= 0; i<ntlfactors.length(); ++i) {
            factors[i] = new GivPolZZDense::Element( deg(ntlfactors[i].a)+1 );
            for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
                NTL::GetCoeff( factors[i]->operator[](j),ntlfactors[i].a,j);
            }
            exp[i] = ntlfactors[i].b;
        }
        return factors;
    } 

#endif

    typedef GivPolynomialRing<Modular<double>,Dense> GivPolMdDense;
    template <>
    template <>
    std::vector<GivPolMdDense::Element *>& 
    GivPolMdDense::factor (std::vector<GivPolMdDense::Element* > & factors, 
                           std::vector<unsigned long>& exp,
                           const GivPolMdDense::Element& P)
    {
	integer charac;
	_domain.characteristic(charac);
	double p = charac;
        typedef GivaroField<Modular<double> > GivModDouble;
        typedef Poly1FactorDom< GivModDouble,Dense, GivModDouble::RandIter> PolysDouble;
        

	PolysDouble PFD(*this, GivModDouble::RandIter(_domain));
	std::vector<PolysDouble::Element> factors2;
	PFD.CZfactor ( factors2, exp, static_cast<PolysDouble::Element>(P),p);

            //std::cerr<<"factorization done"<<std::endl;
	factors.resize(factors2.size());
	std::vector<GivPolMdDense::Element* >::iterator itf = factors.begin();
	std::vector<PolysDouble::Element >::const_iterator itf2 = factors2.begin();
	for (; itf2 != factors2.end();++itf,++itf2){
            *itf = new GivPolMdDense::Element(*itf2);
		//std::cerr<<"converting factor"<<(*itf)<<std::endl;
            for (size_t i=0; i< (*itf)->size();++i)
                _domain.divin((*itf)->operator[](i),(*itf)->operator[]((*itf)->size()-1));
            _domain.assign((*itf)->operator[]((*itf)->size()-1),1.0);
	}
	return factors;
    }	

} // namespace LinBox


#endif // __LINBOX_givaropolynomial_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
