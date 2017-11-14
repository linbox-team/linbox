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
#include "linbox/integer.h"

#include <givaro/givpoly1factor.h>

    // Namespace in which all LinBox code resides
namespace LinBox {
    template<class BaseRing>
    class DensePolynomial;
    
        /** Polynomials.
         * \ingroup ring
         *
         *  @tparam Domain
         *  @tparam StorageTag
         */
	template <class BaseRing, class StorageTag= Givaro::Dense>
	class PolynomialRing : public Givaro::Poly1FactorDom<BaseRing,StorageTag> {
	public:
        
		typedef typename Givaro::Poly1FactorDom<BaseRing,StorageTag> Parent_t;
		typedef DensePolynomial<BaseRing> Element;
		typedef Element Polynomial;
		typedef Element Rep;

        PolynomialRing (const BaseRing& R) : Parent_t(R) {}
        
        PolynomialRing (const BaseRing& R, const Givaro::Indeter& I) : Parent_t(R, I) {}

                   // -- Init polynomial adds his base field
        template<typename... Args>
        Rep& init(Rep& p, Args... args) const {
            Parent_t::init(static_cast<typename Parent_t::Element&>(p),args...);
            p._field = &Parent_t::subdomain();
            return p;
        }

		template<template<class,class> class Vector,template <class> class Alloc>
        Vector<Polynomial, Alloc<Polynomial> >&
        factor (Vector<Polynomial,Alloc<Polynomial> >& factors,
                std::vector<uint64_t>& exp,
                const Polynomial& P){

            Vector<typename Parent_t::Element,Alloc<typename Parent_t::Element> > giv_factors;
            this->CZfactor(giv_factors, exp, P); // Cantor-Zassenhaus factorization
            factors.clear();
            for (size_t i=0; i<giv_factors.size();i++){
                factors.push_back(Element(giv_factors[i],this->_domain));
            }
            return factors;
        }
    };
}

#include "linbox/polynomial/dense-polynomial.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/ring/ntl.h"
#include "NTL/ZZXFactoring.h"

namespace LinBox{
    template<>
    template<>
    std::vector<LinBox::DensePolynomial<Givaro::ZRing<Givaro::Integer> > >&
    LinBox::PolynomialRing<Givaro::ZRing<Givaro::Integer>,Givaro::Dense>::
    factor<std::vector> (std::vector<LinBox::DensePolynomial<Givaro::ZRing<Givaro::Integer> > >& factors,
                         std::vector<uint64_t>& exp,
                         const LinBox::DensePolynomial<Givaro::ZRing<Givaro::Integer> >&P)
    {
        NTL::ZZXFac_InitNumPrimes = 1;
        NTL::ZZX f;
        for (size_t i = 0; i < P.size(); ++i){
            NTL::SetCoeff (f, (long)i, NTL::to_ZZ((std::string( P[(size_t)i] )).c_str()) );
        }
        NTL::vec_pair_ZZX_long ntlfactors;
        NTL::ZZ c;
        NTL::factor (c, ntlfactors, f);

        LinBox::NTL_ZZ NTLIntDom;
        Givaro::ZRing<Givaro::Integer> GivZZ;
        factors.clear();
        exp.resize((size_t)ntlfactors.length());

        for (size_t i= 0; i< (size_t)ntlfactors.length(); ++i) {
            NTL::ZZ t;
            DensePolynomial<Givaro::ZRing<Givaro::Integer> > f (GivZZ, static_cast<size_t>(NTL::deg (ntlfactors [i].a) + 1));
            for(size_t j = 0; j <= (size_t)deg (ntlfactors [i].a); ++j) {
                NTL::GetCoeff (t, ntlfactors [i].a, j);
                NTLIntDom.convert (f[j], t);
            }
            factors.push_back(f);
            exp [i] = ntlfactors[i].b;
        }
        return factors;
    }
} // namespace LinBox


// #include "givaro/zring.h"
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
namespace LinBox{
    template <>
    template <>
    std::vector<DensePolynomial<NTL_ZZ> >&
    PolynomialRing<NTL_ZZ,Givaro::Dense>::factor (std::vector<DensePolynomial<NTL_ZZ> >& factors,
                                                  std::vector<uint64_t>& exp,
                                                  const DensePolynomial<NTL_ZZ> &P)
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
        NTL_ZZ NTLIntDom;
        factors.clear();//resize (ntlfactors.length());
        exp.resize((size_t)ntlfactors.length());
        for (size_t i= 0; i < (size_t)ntlfactors.length(); ++i) {
            DensePolynomial<NTL_ZZ> f (NTLIntDom, static_cast<size_t>(NTL::deg(ntlfactors[i].a)+1));
            for(size_t j = 0; j <= (size_t) deg(ntlfactors[i].a); ++j) {
                NTL::GetCoeff (f[j], ntlfactors [i].a, j);
            }
            factors.push_back(f);
            exp [i] = ntlfactors[i].b;
        }
        return factors;
    }
}

#endif

    // 	typedef GivaroPolynomialRing<Givaro::Modular<double>, Givaro::Dense> GivaroPolMdDense;

// 	template <>
// 	template <>
// 	std::vector<GivaroPolMdDense::Element *>&
// 	GivaroPolMdDense::factor (std::vector<GivaroPolMdDense::Element* > & factors,
// 			       std::vector<uint64_t>& exp,
// 			       const GivaroPolMdDense::Element& P)
// 	{
// 		integer charac;
// 		_domain.characteristic(charac);
// 		double p = charac;
// 		typedef Givaro::Modular<double> GivModDouble;
// 		typedef Givaro::Poly1FactorDom< GivModDouble, Givaro::Dense, GivModDouble::RandIter> PolysDouble;


// 		PolysDouble PFD(*this, GivModDouble::RandIter(_domain));
// 		std::vector<PolysDouble::Element> factors2;
// 		PFD.CZfactor ( factors2, exp, static_cast<PolysDouble::Element>(P),p);

// 		//std::cerr<<"factorization done"<<std::endl;
// 		factors.resize((size_t)factors2.size());
// 		std::vector<GivaroPolMdDense::Element* >::iterator itf = factors.begin();
// 		std::vector<PolysDouble::Element >::const_iterator itf2 = factors2.begin();
// 		for (; itf2 != factors2.end();++itf,++itf2){
// 			*itf = new GivaroPolMdDense::Element(*itf2);
// 			//std::cerr<<"converting factor"<<(*itf)<<std::endl;
// 			for (size_t i=0; i< (*itf)->size();++i)
// 				_domain.divin((*itf)->operator[](i),(*itf)->operator[]((*itf)->size()-1));
// 			_domain.assign((*itf)->operator[]((*itf)->size()-1),1.0);
// 		}
// 		return factors;
// 	}

//} // namespace Givaro


// Dense univariate polynomials are manipulated by Givaro's polynomial factorization domain:
// namespace LinBox{
//     template <class Domain,
//               class Tag = Givaro::Dense,
//               class RandomIterator = Givaro::GivRandom>
//     using PolynomialRing = Givaro::Poly1FactorDom <Domain, Tag, RandomIterator>;
// }

#endif // __LINBOX_givaropolynomial_H


// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

