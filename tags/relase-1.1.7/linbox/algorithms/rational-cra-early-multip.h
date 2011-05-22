/* Copyright (C) 2010 LinBox
 * Written by <Jean-Guillaume.Dumas@imag.fr> 
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

#ifndef __LINBOX_rational_early_multip_cra_H
#define __LINBOX_rational_early_multip_cra_H

#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/rational-cra-early-single.h"
#include "linbox/algorithms/rational-cra-full-multip.h"

namespace LinBox {
	
	template<class Domain_Type>
	struct EarlyMultipRatCRA : public EarlySingleRatCRA<Domain_Type>, public FullMultipRatCRA<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain_Type::Element 	DomainElement;
		typedef EarlyMultipRatCRA<Domain>	Self_t;
	protected:
		// Random coefficients for a linear combination 
		// of the elements to be reconstructed
		std::vector< unsigned long >      	randv;
		
	public:
		
		
		EarlyMultipRatCRA(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD) 
			: EarlySingleRatCRA<Domain>(EARLY), FullMultipRatCRA<Domain>() {
		}
		
		
		template<template<class, class> class Vect, template <class> class Alloc>
		void initialize (const Domain& D, const Vect <DomainElement, Alloc<DomainElement> >& e) {
			// Random coefficients for a linear combination 
			// of the elements to be reconstructed
			srand48(BaseTimer::seed());
			randv. resize ( e.size() );
			for ( std::vector<unsigned long>::iterator int_p = randv. begin(); 
			      int_p != randv. end(); ++ int_p) 
				*int_p = ((unsigned long)lrand48()) % 20000;        
			
			DomainElement z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			EarlySingleRatCRA<Domain>::initialize(D,dot(z, D, e, randv) );
			FullMultipRatCRA<Domain>::initialize(D, e);
		}
		
		
		template<template<class,class> class Vect, template <class> class Alloc> 
		void progress (const Domain& D, const Vect<DomainElement, Alloc<DomainElement> >& e) {
			DomainElement z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			EarlySingleRatCRA<Domain>::progress(D, dot(z, D, e, randv));
			FullMultipRatCRA<Domain>::progress(D, e);
		}
		
		template<template<class, class> class Vect, template <class> class Alloc>
		Vect<Integer, Alloc<Integer> >& result(Vect<Integer, Alloc<Integer> >& num, Integer& den) {
			return FullMultipRatCRA<Domain>::result(num, den);
		}
		
		bool terminated() {
			return EarlySingleRatCRA<Domain>::terminated();
		}
		
		bool noncoprime(const Integer& i) const {
			return EarlySingleRatCRA<Domain>::noncoprime(i);
		}
		
	protected:
		
		template <template<class, class> class Vect1, template <class> class Alloc, class Vect2>
		DomainElement& dot (DomainElement& z, const Domain& D, 
				    const Vect1<DomainElement, Alloc<DomainElement> >& v1, const Vect2& v2) {
			
			D.init(z,0); DomainElement tmp;
			typename Vect1<DomainElement, Alloc<DomainElement> >::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
				for (v1_p  = v1. begin(), v2_p = v2. begin(); 
				     v1_p != v1. end(); 
				     ++ v1_p, ++ v2_p)       
					D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));
            
				//             commentator.report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "v: " << v2 << std::endl;
				//             commentator.report(Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION) << "z: " << z << std::endl;
				return z;
			}
			};
	}
		
#endif //__LINBOX_rational_early_multip_cra_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
