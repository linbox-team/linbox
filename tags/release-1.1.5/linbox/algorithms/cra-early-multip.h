/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Time-stamp: <12 Mar 07 19:38:54 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_CRA_EARLY_MULTIP_H
#define __LINBOX_CRA_EARLY_MULTIP_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

#include "linbox/algorithms/cra-early-single.h"
#include "linbox/algorithms/cra-full-multip.h"



namespace LinBox {
    
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
      
        

	protected:
        
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
#endif

