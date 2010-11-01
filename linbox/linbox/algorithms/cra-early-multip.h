/* Copyright (C) 2007 LinBox
 * Written by JG Dumas
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


#ifndef __LINBOX_cra_early_multip_H
#define __LINBOX_cra_early_multip_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

#include "linbox/algorithms/cra-early-single.h"
#include "linbox/algorithms/cra-full-multip.h"


namespace LinBox 
{
    
	template<class Domain_Type>
	struct EarlyMultipCRA : public EarlySingleCRA<Domain_Type>, public FullMultipCRA<Domain_Type> {
		typedef Domain_Type			Domain;
		typedef typename Domain::Element DomainElement;
		typedef EarlyMultipCRA<Domain> 		Self_t;

	protected:
		// Random coefficients for a linear combination 
		// of the elements to be reconstructed
		std::vector< unsigned long >      	randv;

	public:
	
		EarlyMultipCRA(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD)
			: EarlySingleCRA<Domain>(EARLY), FullMultipCRA<Domain>() {}

        Integer& getModulus(Integer& m) {EarlySingleCRA<Domain>::getModulus(m); return m;}
	Integer& getResidue(Integer& m) {EarlySingleCRA<Domain>::getResidue(m); return m;}
        template<template<class T> class Vect>
        Vect<Integer>& getResidue(Vect<Integer>& m) {FullMultipCRA<Domain>::getResidue(m); return m;}
     
	template<template<class T> class Vect>
        void initialize (const Integer& D, const Vect<Integer>& e) {
            srand48(BaseTimer::seed());
            randv. resize ( e.size() );
            for ( std::vector<unsigned long>::iterator int_p = randv. begin(); int_p != randv. end(); ++ int_p)
	            *int_p = ((unsigned long)lrand48()) % 20000;
            Integer z;
	    dot(z, D, e, randv);
            EarlySingleCRA<Domain>::initialize(D, z);
            FullMultipCRA<Domain>::initialize(D, e);
        }
	
	template<template <class> class Alloc, template<class, class> class Vect>
	void initialize (const Domain& D, const Vect<DomainElement, Alloc<DomainElement> >& e) {
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
		EarlySingleCRA<Domain>::initialize(D,dot(z, D, e, randv) );
		FullMultipCRA<Domain>::initialize(D, e);
	}

	template<template<class T> class Vect>
        void progress (const Integer& D, const Vect<Integer>& e) {

		Integer z;
	        EarlySingleCRA<Domain>::progress(D, dot(z, D, e, randv));
	        FullMultipCRA<Domain>::progress(D, e);
	}

		template<template <class> class Alloc, template<class, class> class Vect>
		void progress (const Domain& D, const Vect<DomainElement, Alloc<DomainElement> >& e) {
			DomainElement z;
			// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			EarlySingleCRA<Domain>::progress(D, dot(z, D, e, randv));
			FullMultipCRA<Domain>::progress(D, e);
		}

		template<template <class> class Alloc, template<class, class> class Vect>
		Vect<Integer, Alloc<Integer> >& result(Vect<Integer, Alloc<Integer> >& d) {
			return FullMultipCRA<Domain>::result(d);
		}

		bool terminated() {
			return EarlySingleCRA<Domain>::terminated();
		}
        
		bool noncoprime(const Integer& i) const {
			return EarlySingleCRA<Domain>::noncoprime(i);
		}
      
        bool changeVector() {
                for ( std::vector<unsigned long>::iterator int_p = randv. begin();int_p != randv. end(); ++ int_p)
	                *int_p = ((unsigned long)lrand48()) % 20000;

		std::vector<Integer> e(randv.size());
		/* clear CRAEarlySingle; */				
		EarlySingleCRA<Domain>::occurency_ = 0;
		EarlySingleCRA<Domain>::nextM_ = 1UL;
		EarlySingleCRA<Domain>::primeProd_ = 1UL;
		EarlySingleCRA<Domain>::residue_ = 0;
	        
		/* Computation of residue_ */
		std::vector< LazyProduct >::iterator _mod_it = FullMultipCRA<Domain>::RadixPrimeProd_.begin();// list of prime products
		std::vector< std::vector<Integer> >::iterator _tab_it = FullMultipCRA<Domain>::RadixResidues_.begin();// list of residues as vectors of size 1
		std::vector< bool >::iterator    _occ_it = FullMultipCRA<Domain>::RadixOccupancy_.begin();//flags of occupied fields
		int prev_shelf=0, shelf = 0;
		for (;_occ_it != FullMultipCRA<Domain>::RadixOccupancy_.end(); ++_mod_it, ++_tab_it, ++_occ_it ) {
			++shelf;
			if (*_occ_it) {
				Integer D = _mod_it->operator()();
				std::vector<Integer> e(randv.size());
				e = *_tab_it;
				Integer z;
				dot(z,D, e, randv);
				Integer prev_residue_ = EarlySingleCRA<Domain>::residue_;
				EarlySingleCRA<Domain>::progress(D,z);
                		if (prev_residue_ == EarlySingleCRA<Domain>::residue_ )
					EarlySingleCRA<Domain>::occurency_ = EarlySingleCRA<Domain>::occurency_ +  (shelf - prev_shelf);
	                        if ( EarlySingleCRA<Domain>::terminated() ) {
					return true;
				}
				prev_shelf = shelf;
			}
		}
		return false;
	}

	protected:
        
		template <template<class T> class Vect1, class Vect2>
	        Integer& dot (Integer& z, const Integer& D, const Vect1<Integer>& v1, const Vect2& v2) {
			z = 0;
			typename Vect1<Integer>::const_iterator v1_p;
			typename Vect2::const_iterator v2_p;
			for (v1_p  = v1. begin(), v2_p = v2. begin(); v1_p != v1. end(); ++ v1_p, ++ v2_p) {
				z = (z + (*v1_p)*(*v2_p))%D;
			}
			return z;
		}

		template <template <class> class Alloc, template<class, class> class Vect1, class Vect2>
		DomainElement& dot (DomainElement& z, const Domain& D, 
				    const Vect1<DomainElement, Alloc<DomainElement> >& v1,
				    const Vect2& v2) {
            
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
#endif //__LINBOX_cra_early_multip_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
