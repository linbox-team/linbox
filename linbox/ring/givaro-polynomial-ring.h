/* linbox/ring/givaro-polynomial.h
 * Copyright(C) LinBox
 * Written by
 * Clement Pernet
 *
 *
 * ========LICENCE========
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
 *.
 */

/*! @file ring/givaro-polynomial-ring.h
 * @ingroup ring
 * @brief NO DOC
 */

#ifndef __LINBOX_givaro_polynomial_ring_H
#define __LINBOX_givaro_polynomial_ring_H

#include <iostream>
#include <givaro/zring.h>
#include "linbox/polynomial/givaro-polynomial.h"
#include "linbox/integer.h"


// Dense univariate polynomials are manipulated by Givaro's polynomial factorization domain:
#include <givaro/givpoly1factor.h>

// This file the redefinitions for some methods:  e.g. NTL is faster for integer polynomial factorization

// // Namespace in which all LinBox code resides
// namespace LinBox
// {
    
    //     /** Polynomials.
    //      * \ingroup ring
    //      *
    //      *  @tparam Domain
    //      *  @tparam StorageTag
    //      */
	// template <class Domain, class StorageTag= Givaro::Dense>
	// class GivaroPolynomialRing : public Givaro::Poly1FactorDom< Domain,StorageTag> {
	// public:
        
	// 	typedef typename Givaro::Poly1FactorDom<Domain,StorageTag> Parent_t;
	// 	typedef GivaroPolynomial<Domain,StorageTag> Element;
	// 	typedef Element Polynomial;

	// 	GivaroPolynomialRing (const Domain& D) : Parent_t( Domain(D) )
	// 	{}
        
	// 	GivaroPolynomialRing (const Domain& D, const Givaro::Indeter& I) :
    //             Parent_t(Domain(D), I)
    //         {}
        
	// 	template< template<class, class> class Container, template <class> class Alloc>
    //     Container<Polynomial,Alloc<Polynomial> >&
    //     factor (Container<Polynomial,Alloc<Polynomial> >& factors, std::vector<uint64_t>& exp, const Polynomial& P){
    //         CZfactor(factors, exp, P); // Cantor-Zassenhaus factorization
    //         return factors;
    //     }
        

	// };


#ifdef __LINBOX_HAVE_NTL
//}

#include "linbox/ring/ntl.h"
#include "NTL/ZZXFactoring.h"
namespace LinBox
{
	typedef GivaroPolynomialRing<Givaro::ZRing<integer>, Givaro::Dense> GivaroPolIntDense;

	template <>
	template <>
	std::vector<GivaroPolynomial<GivaroPolIntDense::Element> >&
	GivaroPolIntDense::factor (std::vector<GivaroPolIntDense::Element>& factors,
                               std::vector<uint64_t>& exp,
                               const  &P)
	{
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < P.size(); ++i){
			NTL::SetCoeff (f, (long)i, NTL::to_ZZ((std::string( P[(size_t)i] )).c_str()) );
		}
		NTL::vec_pair_ZZX_long ntlfactors;
		NTL::ZZ c;
		NTL::factor (c, ntlfactors, f);

		NTL::ZZ t;
		NTL_ZZ NTLIntDom;
		factors.resize((size_t)ntlfactors.length());
		exp.resize((size_t)ntlfactors.length());
		for (int i= 0; i<ntlfactors.length(); ++i) {
			factors[(size_t)i] = new GivaroPolIntDense::Element( (size_t)( deg(ntlfactors[i].a)+1 ));
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff(t,ntlfactors[i].a,j);
				NTLIntDom.convert( factors[(size_t)i]->operator[]((size_t)j), t );
			}
			exp[(size_t)i] =(unsigned long) ntlfactors[i].b;
		}
		return factors;
	}
}

// #include "givaro/zring.h"
 namespace LinBox
 {
// 	typedef GivaroPolynomialRing<Givaro::ZRing<Integer>, Givaro::Dense> GivaroPolPIDIntDense;
// 	template <>
// 	template <>
// 	std::vector<GivaroPolPIDIntDense::Element* >&
// 	GivaroPolPIDIntDense::factor<std::vector<GivaroPolPIDIntDense::Element* > > (std::vector<GivaroPolPIDIntDense::Element* >& factors,
// 			std::vector<unsigned long>& exp,
// 			const GivaroPolPIDIntDense::Element &P)
// 	{
// 		NTL::ZZXFac_InitNumPrimes = 1;
// 		NTL::ZZX f;
// 		for (size_t i = 0; i < P.size(); ++i){
// 			NTL::SetCoeff (f, (long)i, NTL::to_ZZ((std::string( P[(size_t)i] )).c_str()) );
// 		}
// 		NTL::vec_pair_ZZX_long ntlfactors;
// 		NTL::ZZ c;
// 		NTL::factor (c, ntlfactors, f);

// 		NTL::ZZ t;
// 		NTL_ZZ NTLIntDom;
// 		factors.resize((size_t)ntlfactors.length());
// 		exp.resize((size_t)ntlfactors.length());
// 		for (int i= 0; i<ntlfactors.length(); ++i) {
// 			factors[(size_t)i] = new GivaroPolPIDIntDense::Element( (size_t)deg(ntlfactors[i].a)+1 );
// 			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
// 				NTL::GetCoeff(t,ntlfactors[i].a,j);
// 				NTLIntDom.convert( factors[(size_t)i]->operator[]((size_t)j), t );
// 			}
// 			exp[(size_t)i] = (unsigned long)ntlfactors[i].b;
// 		}
// 		return factors;
// 	}

	typedef GivaroPolynomialRing< NTL_ZZ , Givaro::Dense> GivaroPolZZDense;
	template <>
	template <>
	std::vector<GivaroPolZZDense::Element* >&
	GivaroPolZZDense::factor<std::vector<GivaroPolZZDense::Element* > >
	(std::vector<GivaroPolZZDense::Element* >& factors,
	 std::vector<uint64_t>& exp,
	 const GivaroPolZZDense::Element &P)
	{
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < P.size(); ++i){
			NTL::SetCoeff (f, (long)i, P[(size_t)i]);
		}
		NTL::vec_pair_ZZX_long ntlfactors;
		NTL::ZZ c;
		NTL::factor (c, ntlfactors, f);

		NTL::ZZ t;
		factors.resize((size_t)ntlfactors.length());
		exp.resize((size_t)ntlfactors.length());
		for (int i= 0; i<ntlfactors.length(); ++i) {
			factors[(size_t)i] = new GivaroPolZZDense::Element((size_t) deg(ntlfactors[i].a)+1 );
			for(int j = 0; j <= deg(ntlfactors[i].a); ++j) {
				NTL::GetCoeff( factors[(size_t)i]->operator[]((size_t)j),ntlfactors[i].a,j);
			}
			exp[(size_t)i] = (unsigned long)ntlfactors[i].b;
		}
		return factors;
	}

#endif

	typedef GivaroPolynomialRing<Givaro::Modular<double>, Givaro::Dense> GivaroPolMdDense;

	template <>
	template <>
	std::vector<GivaroPolMdDense::Element *>&
	GivaroPolMdDense::factor (std::vector<GivaroPolMdDense::Element* > & factors,
			       std::vector<uint64_t>& exp,
			       const GivaroPolMdDense::Element& P)
	{
		integer charac;
		_domain.characteristic(charac);
		double p = charac;
		typedef Givaro::Modular<double> GivModDouble;
		typedef Givaro::Poly1FactorDom< GivModDouble, Givaro::Dense, GivModDouble::RandIter> PolysDouble;


		PolysDouble PFD(*this, GivModDouble::RandIter(_domain));
		std::vector<PolysDouble::Element> factors2;
		PFD.CZfactor ( factors2, exp, static_cast<PolysDouble::Element>(P),p);

		//std::cerr<<"factorization done"<<std::endl;
		factors.resize((size_t)factors2.size());
		std::vector<GivaroPolMdDense::Element* >::iterator itf = factors.begin();
		std::vector<PolysDouble::Element >::const_iterator itf2 = factors2.begin();
		for (; itf2 != factors2.end();++itf,++itf2){
			*itf = new GivaroPolMdDense::Element(*itf2);
			//std::cerr<<"converting factor"<<(*itf)<<std::endl;
			for (size_t i=0; i< (*itf)->size();++i)
				_domain.divin((*itf)->operator[](i),(*itf)->operator[]((*itf)->size()-1));
			_domain.assign((*itf)->operator[]((*itf)->size()-1),1.0);
		}
		return factors;
	}

} // namespace LinBox


#endif // __LINBOX_givaropolynomial_H


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

