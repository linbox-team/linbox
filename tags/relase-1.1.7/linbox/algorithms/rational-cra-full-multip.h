/* Copyright (C) 2007  LinBox
 * Written by JG Dumas
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_rational_full_multip_cra_H
#define __LINBOX_rational_full_multip_cra_H

#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/cra-full-multip.h"

namespace LinBox 
{

#if 0
	template<class T, template <class T> class Container>
		std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
			for(typename Container<T>::const_iterator refs =  C.begin();
					refs != C.end() ;
					++refs )
				o << (*refs) << " " ;
			return o << std::endl;
		}
#endif

    template<class Domain_Type>
    struct FullMultipRatCRA : public virtual FullMultipCRA<Domain_Type> {
        typedef Domain_Type				Domain;
        typedef FullMultipCRA<Domain> 			Father_t;
        typedef typename Father_t::DomainElement 	DomainElement;
        typedef FullMultipRatCRA<Domain>		Self_t;
        PID_integer _ZZ;
    public:

        using Father_t::RadixSizes_;
        using Father_t::RadixResidues_; 
        using Father_t::RadixPrimeProd_;
        using Father_t::RadixOccupancy_;
        

        FullMultipRatCRA(const double BOUND = 0.0) : Father_t(BOUND) {}

		
	    template<template<class, class> class Vect, template <class> class Alloc>
	    Vect<Integer, Alloc<Integer> >& result (Vect<Integer, Alloc<Integer> > &num, Integer& den){
            num.resize( (Father_t::RadixResidues_.front()).size() );
            std::vector< LazyProduct >::iterator 			_mod_it = Father_t::RadixPrimeProd_.begin();
            std::vector< std::vector< Integer > >::iterator _tab_it = Father_t::RadixResidues_.begin();
            std::vector< bool >::iterator    				_occ_it = Father_t::RadixOccupancy_.begin();
            LazyProduct Product;
            for( ; _occ_it != Father_t::RadixOccupancy_.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    Product = *_mod_it;
                    std::vector<Integer>::iterator t0_it = num.begin();
                    std::vector<Integer>::iterator t_it = _tab_it->begin();
                    if (++_occ_it == Father_t::RadixOccupancy_.end()) {
                        den = 1;
                        Integer s, nd; _ZZ.sqrt(s, _mod_it->operator()());
                        for( ; t0_it != num.end(); ++t0_it, ++t_it) {
                            iterativeratrecon(*t0_it = *t_it, nd, den, _mod_it->operator()(), s);
                            if (nd > 1) {
                                std::vector<Integer>::iterator  t02 = num.begin();
                                for( ; t02 != t0_it ; ++t02)
                                    *t02 *= nd;
                                den *= nd;
                            }
                        }
                        return num;
                    } else {
                        for( ; t0_it != num.end(); ++t0_it, ++t_it)
                            *t0_it  = *t_it;
                        ++_mod_it; ++_tab_it; 
                        break;
                    }
                }
            }
            for( ; _occ_it != Father_t::RadixOccupancy_.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    std::vector<Integer>::iterator t0_it = num.begin();
                    std::vector<Integer>::const_iterator t_it = _tab_it->begin();
		    Integer invprod;
                    this->precomputeInvProd(invprod, Product(), _mod_it->operator()() );
                    for( ; t0_it != num.end(); ++t0_it, ++t_it)
                        this->smallbigreconstruct(*t0_it, *t_it, invprod );
                    Product.mulin(*_mod_it);

		    // Moding out and normalization
                    for(t0_it = num.begin();t0_it != num.end(); ++t0_it) {
                        *t0_it %= Product();
                        Integer tmp(*t0_it);
                        this->normalize(*t0_it, tmp, Product());
                    }
                }
            }
            den = 1;
            Integer s, nd; _ZZ.sqrt(s, Product.operator()());
            std::vector<Integer>::iterator t0_it = num.begin();
            for( ; t0_it != num.end(); ++t0_it) {
                iterativeratrecon(*t0_it, nd, den, Product.operator()(), s);
                if (nd > 1) {
                    std::vector<Integer>::iterator  t02 = num.begin();
                    for( ; t02 != t0_it ; ++t02)
                        *t02 *= nd;
                    den *= nd;
                }
            }
            return num;
        }
		

    protected:
        Integer& iterativeratrecon(Integer& u1, Integer& new_den, const Integer& old_den, const Integer& m1, const Integer& s) {
            Integer a;
            _ZZ.reconstructRational(a, new_den, u1*=old_den, m1, s);
            return u1=a;
        }
    };
}

#endif //__LINBOX_rational_full_multip_cra_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
